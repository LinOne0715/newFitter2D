#include "cdc2dfit_kernel.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

struct TrackKey {
    int evt = -1;
    int trkId = -1;
    bool operator<(const TrackKey& other) const {
        return std::tie(evt, trkId) < std::tie(other.evt, other.trkId);
    }
};

struct InputHit {
    bool valid = false;
    int evt = -1;
    int trkId = -1;
    int sl = -1;
    int ax = -1;
    int lay = -1;
    int wire = -1;
    int lr = 0;
    int priorityTime = 0;
    int T0 = 9999;
    int driftTick = 0;
    int finderCharge = 0;
};

struct TrackInfo {
    std::array<InputHit, 5> hits{};
    std::array<bool, 5> hasHit{{false, false, false, false, false}};
    int finderCharge = 0;
    int T0 = 9999;
};

static inline std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static inline std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) out.push_back(trim(tok));
    return out;
}

static inline std::map<std::string, int> build_header_index(const std::vector<std::string>& header) {
    std::map<std::string, int> idx;
    for (int i = 0; i < (int)header.size(); ++i) idx[trim(header[i])] = i;
    return idx;
}

static inline std::string get_field(
    const std::vector<std::string>& fields,
    const std::map<std::string, int>& idx,
    const std::string& key,
    const std::string& def = "")
{
    auto it = idx.find(key);
    if (it == idx.end()) return def;
    int pos = it->second;
    if (pos < 0 || pos >= (int)fields.size()) return def;
    return fields[pos];
}

static inline int stoi_safe(const std::string& s, int def = 0) {
    try { return std::stoi(s); } catch (...) { return def; }
}

static inline int sl_to_ax_index(int sl) {
    switch (sl) {
        case 0: return 0;
        case 2: return 1;
        case 4: return 2;
        case 6: return 3;
        case 8: return 4;
        default: return -1;
    }
}

static inline bool load_hits_with_t0_csv(
    const std::string& path,
    std::map<TrackKey, TrackInfo>& tracks)
{
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "[csv] cannot open " << path << "\n";
        return false;
    }

    std::string line;
    if (!std::getline(fin, line)) {
        std::cerr << "[csv] empty file " << path << "\n";
        return false;
    }

    auto header = split_csv_line(line);
    auto idx = build_header_index(header);

    const std::vector<std::string> required = {
        "evt", "trkId", "sl", "lay", "wire", "lr",
        "priorityTime", "T0", "driftTick", "finderCharge"
    };
    for (const auto& key : required) {
        if (!idx.count(key)) {
            std::cerr << "[csv] missing required column: " << key << "\n";
            return false;
        }
    }

    while (std::getline(fin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto fields = split_csv_line(line);

        InputHit h;
        h.evt          = stoi_safe(get_field(fields, idx, "evt", "-1"), -1);
        h.trkId        = stoi_safe(get_field(fields, idx, "trkId", "-1"), -1);
        h.sl           = stoi_safe(get_field(fields, idx, "sl", "-1"), -1);
        h.ax           = sl_to_ax_index(h.sl);
        h.lay          = stoi_safe(get_field(fields, idx, "lay", "-1"), -1);
        h.wire         = stoi_safe(get_field(fields, idx, "wire", "-1"), -1);
        h.lr           = stoi_safe(get_field(fields, idx, "lr", "0"), 0);
        h.priorityTime = stoi_safe(get_field(fields, idx, "priorityTime", "0"), 0);
        h.T0           = stoi_safe(get_field(fields, idx, "T0", "9999"), 9999);
        h.driftTick    = stoi_safe(get_field(fields, idx, "driftTick", "0"), 0);
        h.finderCharge = stoi_safe(get_field(fields, idx, "finderCharge", "0"), 0);
        h.valid        = true;

        if (h.evt < 0 || h.trkId < 0 || h.ax < 0) continue;

        auto& info = tracks[{h.evt, h.trkId}];
        info.hasHit[h.ax] = true;
        info.hits[h.ax] = h;
        info.finderCharge = h.finderCharge;
        info.T0 = h.T0;
    }

    return true;
}

int main(int argc, char** argv)
{
    const std::string hitsPath = (argc > 1) ? argv[1] : "hits_with_t0.csv";
    const std::string outPath  = (argc > 2) ? argv[2] : "hls_local_aligned_out_nts.csv";

    std::map<TrackKey, TrackInfo> tracks;
    if (!load_hits_with_t0_csv(hitsPath, tracks)) return 1;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot write " << outPath << "\n";
        return 1;
    }
    out << std::setprecision(15);
    out << "evt,trkId,phi0,omega,chi2,nhits,rho,chargeFit,phi0_internal,T0\n";

    int printed = 0;
    for (const auto& kv : tracks) {
        const TrackKey& key = kv.first;
        const TrackInfo& info = kv.second;

        ap_uint<6>  layer[5];
        ap_uint<4>  sl_raw[5];
        ap_uint<10> wire[5];
        ap_uint<2>  lr[5];
        ap_uint<10> priorityTime[5];
        ap_uint<14> driftTick[5];

        for (int i = 0; i < 5; ++i) {
            layer[i] = 63;
            sl_raw[i] = 15;
            wire[i] = 0;
            lr[i] = 0;
            priorityTime[i] = 0;
            driftTick[i] = 0;
        }

        for (int i = 0; i < 5; ++i) {
            if (!info.hasHit[i]) continue;
            const auto& h = info.hits[i];
            layer[i] = (unsigned)h.lay;
            sl_raw[i] = (unsigned)h.sl;
            wire[i] = (unsigned)h.wire;
            lr[i] = (unsigned)h.lr;
            priorityTime[i] = (unsigned)std::max(0, h.priorityTime);
            int dt = h.driftTick;
            if (dt < 0) dt = 0;
            if (dt > 511) dt = 511;
            driftTick[i] = (unsigned)dt;
        }

        fixed_type_general rho = 0;
        fixed_type_general phi0_internal = 0;
        fixed_type_general phi0 = 0;
        fixed_type_general omega = 0;
        fixed_type_general chi2 = 0;
        ap_int<3> chargeFit = 0;
        ap_uint<3> nhits = 0;

        cdc2dfit_kernel(
            layer[0], sl_raw[0], wire[0], lr[0], priorityTime[0], driftTick[0],
            layer[1], sl_raw[1], wire[1], lr[1], priorityTime[1], driftTick[1],
            layer[2], sl_raw[2], wire[2], lr[2], priorityTime[2], driftTick[2],
            layer[3], sl_raw[3], wire[3], lr[3], priorityTime[3], driftTick[3],
            layer[4], sl_raw[4], wire[4], lr[4], priorityTime[4], driftTick[4],
            info.T0,
            info.finderCharge,
            rho,
            phi0_internal,
            phi0,
            omega,
            chi2,
            chargeFit,
            nhits
        );

        if (printed < 30) {
            std::cout << std::fixed << std::setprecision(12)
                      << "evt=" << key.evt
                      << " trk=" << key.trkId
                      << " phi0=" << (double)phi0
                      << " omega=" << (double)omega
                      << " chi2=" << (double)chi2
                      << " nhits=" << (unsigned)nhits
                      << " rho=" << (double)rho
                      << " q=" << (int)chargeFit
                      << " T0=" << info.T0
                      << "\n";
            ++printed;
        }

        out << key.evt << ','
            << key.trkId << ','
            << (double)phi0 << ','
            << (double)omega << ','
            << (double)chi2 << ','
            << (unsigned)nhits << ','
            << (double)rho << ','
            << (int)chargeFit << ','
            << (double)phi0_internal << ','
            << info.T0 << '\n';
    }

    std::cout << "Wrote " << tracks.size() << " tracks to " << outPath << "\n";
    return 0;
}
