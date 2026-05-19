#!/usr/bin/env python3
import basf2
import ROOT
import csv
import sys

SL_TO_LAYER = {0: 3, 2: 16, 4: 28, 6: 40, 8: 52}
AXIAL_SL = [0, 2, 4, 6, 8]


class Export2DFitterInputs(basf2.Module):
    def __init__(self, out_csv="hits_with_t0.csv"):
        super().__init__()
        self.out_csv = out_csv
        self.rows = []

    def event(self):
        finder_tracks = ROOT.Belle2.PyStoreArray("TRGCDC2DFinderTracks")
        evt_t0_obj = ROOT.Belle2.PyStoreObj("BinnedEventT0")
        evt_meta = ROOT.Belle2.PyStoreObj("EventMetaData")

        evt = int(evt_meta.obj().getEvent() - 1)

        T0 = 9999
        try:
            if evt_t0_obj.obj().hasBinnedEventT0(ROOT.Belle2.Const.CDC):
                T0 = int(evt_t0_obj.obj().getBinnedEventT0(ROOT.Belle2.Const.CDC))
        except Exception:
            T0 = 9999

        for itrack, trk in enumerate(finder_tracks):
            # same convention as your fitter export
            if hasattr(trk, "getUniqueID"):
                try:
                    trk_id = int(trk.getUniqueID())
                except Exception:
                    trk_id = itrack
            else:
                trk_id = itrack

            finder_charge = 0
            if hasattr(trk, "getChargeSign"):
                try:
                    finder_charge = int(trk.getChargeSign())
                except Exception:
                    finder_charge = 0

            try:
                rel_hits = trk.getRelationsTo("CDCTriggerSegmentHits")
            except Exception:
                continue

            best_hit = {}

            for ihit in range(rel_hits.size()):
                try:
                    weight = rel_hits.weight(ihit)
                except Exception:
                    weight = 1.0

                if weight <= 0:
                    continue

                h = rel_hits[ihit]

                if h.getPriorityPosition() != 3:
                    continue

                iSL = int(h.getISuperLayer())
                if iSL % 2 != 0:
                    continue
                if iSL not in SL_TO_LAYER:
                    continue

                pt = int(h.priorityTime())

                if iSL in best_hit and pt >= best_hit[iSL]["priorityTime"]:
                    continue

                drift_tick = pt - T0
                if drift_tick < 0:
                    drift_tick = 0
                if drift_tick > 511:
                    drift_tick = 511

                best_hit[iSL] = {
                    "evt": evt,
                    "trkId": trk_id,
                    "sl": iSL,
                    "lay": SL_TO_LAYER[iSL],
                    "wire": int(h.getIWire()),
                    "lr": int(h.getLeftRight()),
                    "priorityTime": pt,
                    "T0": T0,
                    "driftTick": drift_tick,
                    "finderCharge": finder_charge,
                }

            for iSL in AXIAL_SL:
                if iSL in best_hit:
                    self.rows.append(best_hit[iSL])

    def terminate(self):
        with open(self.out_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "evt", "trkId", "sl", "lay", "wire", "lr",
                    "priorityTime", "T0", "driftTick", "finderCharge"
                ],
            )
            writer.writeheader()
            for r in self.rows:
                writer.writerow(r)

        print(f"✅ CSV written to {self.out_csv}")
        print(f"Selected axial hits exported: {len(self.rows)}")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: basf2 export_hits_with_t0_module.py <input.root> [out.csv]")

    in_root = sys.argv[1]
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "hits_with_t0.csv"

    main_path = basf2.create_path()
    main_path.add_module("RootInput", inputFileNames=[in_root])
    main_path.add_module(Export2DFitterInputs(out_csv))
    basf2.process(main_path)
    print(basf2.statistics)


if __name__ == "__main__":
    main()