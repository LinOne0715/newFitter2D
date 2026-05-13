#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>

static constexpr std::array<unsigned, 5> PRIORITY_LAYERS = {3, 16, 28, 40, 52};
static constexpr std::array<unsigned, 5> NWIRES          = {160, 192, 256, 320, 384};
static constexpr std::array<double,   5> RADIUS          = {19.8, 40.16, 62.0, 83.84, 105.68};

static constexpr std::array<double, 5> WIRE_PHI_ERR = {
    0.00085106, 0.00039841, 0.00025806, 0.00019084, 0.0001514
};
static constexpr std::array<double, 5> DRIFT_PHI_ERR = {
    0.00085106, 0.00039841, 0.00025806, 0.00019084, 0.0001514
};

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
    int sl = -1;          // raw SL from CSV: 0,2,4,6,8
    int ax = -1;          // mapped index: 0..4
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
    while (std::getline(ss, tok, ',')) {
        out.push_back(trim(tok));
    }
    return out;
}

static inline std::map<std::string, int> build_header_index(const std::vector<std::string>& header) {
    std::map<std::string, int> idx;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        idx[trim(header[i])] = i;
    }
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
    if (pos < 0 || pos >= static_cast<int>(fields.size())) return def;
    return fields[pos];
}

static inline int stoi_safe(const std::string& s, int def = 0) {
    try { return std::stoi(s); } catch (...) { return def; }
}

static inline double stod_safe(const std::string& s, double def = 0.0) {
    try { return std::stod(s); } catch (...) { return def; }
}

static inline double wrap_to_npi_pi(double phi) {
    while (phi > M_PI)  phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
}

static inline double stored_track_phi0(double phi0_internal, double chargeFit) {
    return std::remainder(phi0_internal + chargeFit * M_PI_2, 2.0 * M_PI);
}

static inline double exported_phi0_like_basf2(double phi0_internal) {
    return wrap_to_npi_pi(phi0_internal);
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

static inline bool load_xt_table(
    const std::string& path,
    std::array<std::array<double, 512>, 5>& xt)
{
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "[xt] cannot open " << path << "\n";
        return false;
    }

    std::string line;
    std::size_t row = 0;
    while (std::getline(fin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (row >= 5) {
            std::cerr << "[xt] too many rows in " << path << "\n";
            return false;
        }

        auto toks = split_csv_line(line);
        if (toks.size() != 512) {
            std::cerr << "[xt] row " << row << " has " << toks.size()
                      << " columns, expected 512\n";
            return false;
        }

        for (std::size_t col = 0; col < 512; ++col) {
            xt[row][col] = stod_safe(toks[col], 0.0);
        }
        ++row;
    }

    if (row != 5) {
        std::cerr << "[xt] got " << row << " rows, expected 5\n";
        return false;
    }
    return true;
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
        h.ax           = sl_to_ax_index(h.sl);   // <-- critical fix
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

// ===== basf2-aligned helpers =====

static inline double calPhi_basf2(double wirePhi, double driftLength, double rr, int lr)
{
    double result = wirePhi;
    double driftPhi = std::atan(driftLength / rr);
    if (lr == 1) result -= driftPhi;
    else if (lr == 2) result += driftPhi;
    return result;
}

static inline void rPhiFitter_basf2(
    const std::array<double, 5>& rr,
    const std::array<double, 5>& phi2,
    const std::array<double, 5>& invphierror,
    double& rho,
    double& myphi0,
    double& chi2)
{
    double A = 0.0, B = 0.0, C = 0.0, D = 0.0, E = 0.0;
    std::array<double, 5> invFiterror{};

    for (unsigned i = 0; i < 5; ++i) {
        invFiterror[i] = invphierror[i] / rr[i];
    }

    for (unsigned i = 0; i < 5; ++i) {
        A += std::cos(phi2[i]) * std::cos(phi2[i]) * (invFiterror[i] * invFiterror[i]);
        B += std::sin(phi2[i]) * std::sin(phi2[i]) * (invFiterror[i] * invFiterror[i]);
        C += std::cos(phi2[i]) * std::sin(phi2[i]) * (invFiterror[i] * invFiterror[i]);
        D += rr[i] * std::cos(phi2[i]) * (invFiterror[i] * invFiterror[i]);
        E += rr[i] * std::sin(phi2[i]) * (invFiterror[i] * invFiterror[i]);
    }

    double hcx = D * B - E * C;
    hcx /= 2.0 * (A * B - C * C);

    double hcy = E * A - D * C;
    hcy /= 2.0 * (A * B - C * C);

    rho = std::sqrt(hcx * hcx + hcy * hcy);
    myphi0 = std::atan2(hcy, hcx);
    if (myphi0 < 0.0) myphi0 += 2.0 * M_PI;

    int nTSHits = 0;
    for (unsigned i = 0; i < 5; ++i) {
        if (invphierror[i] != 0.0) ++nTSHits;
    }

    chi2 = 0.0;
    for (unsigned i = 0; i < 5; ++i) {
        if (invphierror[i] == 0.0) continue;
        const double res = 2.0 * (hcx * std::cos(phi2[i]) + hcy * std::sin(phi2[i])) - rr[i];
        chi2 += (res * res) / (invFiterror[i] * invFiterror[i]);
    }

    if (nTSHits > 2) chi2 /= (nTSHits - 2);
    else chi2 = 0.0;
}

static inline void chargeFinder_basf2(
    const std::array<double, 5>& nTSs,
    const std::array<double, 5>& tsIds,
    const std::array<double, 5>& tsHitmap,
    double phi0,
    double inCharge,
    double& foundCharge)
{
    std::array<double, 5> tsPhi{};
    std::array<double, 5> dPhi{};
    std::array<double, 5> dPhi_c{};
    std::array<double, 5> charge{};

    for (unsigned iAx = 0; iAx < 5; ++iAx) {
        tsPhi[iAx] = tsIds[iAx] * 2.0 * M_PI / nTSs[iAx];
        dPhi[iAx] = tsPhi[iAx] - phi0;

        if (dPhi[iAx] < -M_PI) dPhi_c[iAx] = dPhi[iAx] + 2.0 * M_PI;
        else if (dPhi[iAx] > M_PI) dPhi_c[iAx] = dPhi[iAx] - 2.0 * M_PI;
        else dPhi_c[iAx] = dPhi[iAx];

        if (tsHitmap[iAx] == 0.0) charge[iAx] = 0.0;
        else if (dPhi_c[iAx] == 0.0) charge[iAx] = 0.0;
        else if (dPhi_c[iAx] > 0.0) charge[iAx] = 1.0;
        else charge[iAx] = -1.0;
    }

    double sumCharge = charge[0] + charge[1] + charge[2] + charge[3] + charge[4];

    if (sumCharge == 0.0) foundCharge = inCharge;
    else if (sumCharge > 0.0) foundCharge = 1.0;
    else foundCharge = -1.0;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr
            << "Usage: " << argv[0] << " hits_with_t0.csv xt_fixed_5x512.csv\n";
        return 1;
    }

    const std::string hitsPath = argv[1];
    const std::string xtPath   = argv[2];

    std::array<std::array<double, 512>, 5> xt{};
    if (!load_xt_table(xtPath, xt)) return 1;

    std::map<TrackKey, TrackInfo> tracks;
    if (!load_hits_with_t0_csv(hitsPath, tracks)) return 1;

    std::cout << "[initialize] nWires:";
    for (auto w : NWIRES) std::cout << ' ' << w;
    std::cout << "\n[initialize] rr:";
    for (auto r : RADIUS) std::cout << ' ' << r;
    std::cout << "\n\n";

    std::ofstream out("cdc2dfit_out.csv");
    out << "evt,trkId,phi0,omega,chi2,nhits,phi0_internal,phi0_track,chargeFit,T0\n";

    for (const auto& kv : tracks) {
        const TrackKey& key = kv.first;
        const TrackInfo& info = kv.second;

        std::array<double, 5> tsId{};
        std::array<double, 5> wirePhi{};
        std::array<int,    5> priorityTime{};
        std::array<int,    5> driftTick{};
        std::array<double, 5> driftLength{};
        std::array<unsigned, 5> LR{};
        std::array<int, 5> hitIds{};
        std::array<double, 5> useSL{};
        std::array<double, 5> phi2DInvError{};
        std::array<double, 5> phi2D{};

        hitIds.fill(-1);

        unsigned nHits = 0;
        const int T0 = info.T0;
        const double inCharge = static_cast<double>(info.finderCharge);

        for (int i = 0; i < 5; ++i) {
            if (!info.hasHit[i]) continue;

            const auto& h = info.hits[i];
            ++nHits;

            useSL[i] = 1.0;
            hitIds[i] = i;
            tsId[i] = static_cast<double>(h.wire);
            wirePhi[i] = static_cast<double>(h.wire) * 2.0 * M_PI / static_cast<double>(NWIRES[i]);
            LR[i] = static_cast<unsigned>(h.lr);
            priorityTime[i] = h.priorityTime;

            int t = h.driftTick;
            if (t < 0) t = 0;
            if (t > 511) t = 511;
            driftTick[i] = t;

            driftLength[i] = xt[i][t];
        }

        if (nHits < 2) continue;

        for (int i = 0; i < 5; ++i) {
            if (LR[i] == 0) {
                phi2DInvError[i] = 0.0;
                continue;
            }

            if (LR[i] != 3 && T0 != 9999) phi2DInvError[i] = 1.0 / DRIFT_PHI_ERR[i];
            else phi2DInvError[i] = 1.0 / WIRE_PHI_ERR[i];
        }

        for (int i = 0; i < 5; ++i) {
            if (hitIds[i] < 0) continue;

            if (T0 == 9999) phi2D[i] = wirePhi[i];
            else phi2D[i] = calPhi_basf2(wirePhi[i], driftLength[i], RADIUS[i], static_cast<int>(LR[i]));
        }

        std::cout << "=== 2D-fitter debug =========================================\n";
        std::cout << "trkIdx=" << key.trkId << "  nHits=" << nHits << "\n";

        std::cout << "  tsId        :";
        for (double v : tsId) std::cout << ' ' << std::fixed << std::setprecision(6) << v;
        std::cout << "\n";

        std::cout << "  wirePhi(rad):";
        for (double v : wirePhi) std::cout << ' ' << std::fixed << std::setprecision(6) << v;
        std::cout << "\n";

        std::cout << "  priorityTime:";
        for (int v : priorityTime) std::cout << ' ' << v;
        std::cout << "\n";

        std::cout << "  driftTick   :";
        for (int v : driftTick) std::cout << ' ' << v;
        std::cout << "\n";

        std::cout << "  driftLen(cm):";
        for (double v : driftLength) std::cout << ' ' << std::fixed << std::setprecision(6) << v;
        std::cout << "\n";

        std::cout << "  LR          :";
        for (unsigned v : LR) std::cout << ' ' << v;
        std::cout << "\n";

        std::cout << "  hitIds      :";
        for (int v : hitIds) std::cout << ' ' << v;
        std::cout << "\n";

        std::cout << "  useSL(flag) :";
        for (double v : useSL) std::cout << ' ' << v;
        std::cout << "\n";

        std::cout << "  T0          : " << T0 << "\n";
        std::cout << "  finderCharge: " << inCharge << "\n";
        std::cout << "----------------------------------------------------------------\n";

        double rho = 0.0;
        double phi0_internal = 0.0;
        double chi2 = 0.0;

        rPhiFitter_basf2(RADIUS, phi2D, phi2DInvError, rho, phi0_internal, chi2);

        double chargeFit = 0.0;
        std::array<double, 5> nTSs{};
        for (int i = 0; i < 5; ++i) nTSs[i] = static_cast<double>(NWIRES[i]);

        chargeFinder_basf2(nTSs, tsId, useSL, phi0_internal, inCharge, chargeFit);

        double omega = 0.0;
        if (rho != 0.0) omega = chargeFit / rho;

        const double phi0_track = stored_track_phi0(phi0_internal, chargeFit);
        const double phi0_export = exported_phi0_like_basf2(phi0_internal);

        std::cout << std::fixed << std::setprecision(6)
                  << "[fit] rho=" << rho
                  << " phi0_internal=" << phi0_internal
                  << " phi0_track=" << phi0_track
                  << " phi0_export=" << phi0_export
                  << " chi2=" << chi2
                  << " chargeFit=" << chargeFit
                  << " omega=" << omega
                  << "\n\n";

        out << key.evt << ','
            << key.trkId << ','
            << std::setprecision(12) << phi0_export << ','
            << std::setprecision(12) << omega << ','
            << std::setprecision(6)  << chi2 << ','
            << nHits << ','
            << std::setprecision(12) << phi0_internal << ','
            << std::setprecision(12) << phi0_track << ','
            << std::setprecision(1)  << chargeFit << ','
            << T0 << '\n';
    }

    std::cout << "[done] wrote cdc2dfit_out.csv\n";
    return 0;
}