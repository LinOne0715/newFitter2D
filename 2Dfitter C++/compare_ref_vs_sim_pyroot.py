#!/usr/bin/env python3
import os
import math
import argparse
import pandas as pd
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)


def wrap_to_pi(x):
    while x >= math.pi:
        x -= 2.0 * math.pi
    while x < -math.pi:
        x += 2.0 * math.pi
    return x


def rad_to_deg(x):
    return x * 180.0 / math.pi


def choose_text_pos(pos):
    if pos == "upper_right_far":
        return 0.68, 0.88, 0.040
    if pos == "upper_right":
        return 0.62, 0.88, 0.040
    if pos == "upper_left":
        return 0.16, 0.88, 0.040
    return 0.16, 0.88, 0.040


def make_hist_and_graph(values, hist_name, nbins, xmin, xmax):
    hist = ROOT.TH1D(hist_name, "", nbins, xmin, xmax)
    hist.Sumw2()

    for v in values:
        if math.isfinite(v):
            hist.Fill(v)

    xvals, yvals = [], []
    exl, exh, eyl, eyh = [], [], [], []

    for ib in range(1, hist.GetNbinsX() + 1):
        y = hist.GetBinContent(ib)
        if y <= 0:
            continue
        x = hist.GetBinCenter(ib)
        bw = 0.5 * hist.GetBinWidth(ib)

        xvals.append(x)
        yvals.append(y)
        exl.append(bw)
        exh.append(bw)
        err = math.sqrt(y)
        eyl.append(err)
        eyh.append(err)

    graph = ROOT.TGraphAsymmErrors(
        len(xvals),
        array("d", xvals), array("d", yvals),
        array("d", exl), array("d", exh),
        array("d", eyl), array("d", eyh)
    )
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetLineWidth(2)

    return hist, graph


def fit_and_draw(values, title, x_title, out_png, hist_name,
                 hist_range, fit_range, nbins, text_pos):
    xmin, xmax = hist_range
    fmin, fmax = fit_range

    hist, graph = make_hist_and_graph(values, hist_name, nbins, xmin, xmax)

    c = ROOT.TCanvas(f"c_{hist_name}", "", 900, 700)
    c.SetMargin(0.12, 0.04, 0.12, 0.06)

    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle("Entries")
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetMaxDigits(3)

    ymax = max(hist.GetMaximum() * 1.35, 1.0)
    hist.SetMaximum(ymax)
    hist.Draw("HIST")
    graph.Draw("P E SAME")

    fit = ROOT.TF1(f"fit_{hist_name}", "gaus", fmin, fmax)
    fit.SetLineWidth(2)
    fit_result = hist.Fit(fit, "RQ0S")

    mean = hist.GetMean()
    mean_err = hist.GetMeanError()
    std = hist.GetStdDev()
    std_err = hist.GetStdDevError()

    fit_mean = float("nan")
    fit_mean_err = float("nan")
    fit_sigma = float("nan")
    fit_sigma_err = float("nan")

    if fit_result and int(fit_result) == 0:
        fit.Draw("SAME")
        fit_mean = fit.GetParameter(1)
        fit_mean_err = fit.GetParError(1)
        fit_sigma = fit.GetParameter(2)
        fit_sigma_err = fit.GetParError(2)

    x0, y0, dy = choose_text_pos(text_pos)
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(dy)
    latex.DrawLatex(x0, y0,           f"Entries = {hist.GetEntries():.0f}")
    latex.DrawLatex(x0, y0 - 0.06,    f"Hist mean = {mean:.6g} #pm {mean_err:.3g}")
    latex.DrawLatex(x0, y0 - 0.12,    f"Hist std = {std:.6g} #pm {std_err:.3g}")
    if math.isfinite(fit_mean):
        latex.DrawLatex(x0, y0 - 0.18, f"Fit mean = {fit_mean:.6g} #pm {fit_mean_err:.3g}")
        latex.DrawLatex(x0, y0 - 0.24, f"Fit #sigma = {fit_sigma:.6g} #pm {fit_sigma_err:.3g}")

    c.SaveAs(out_png)

    stats = {
        "hist": hist,
        "graph": graph,
        "fit": fit,
        "entries": hist.GetEntries(),
        "hist_mean": mean,
        "hist_mean_err": mean_err,
        "hist_std": std,
        "hist_std_err": std_err,
        "fit_mean": fit_mean,
        "fit_mean_err": fit_mean_err,
        "fit_sigma": fit_sigma,
        "fit_sigma_err": fit_sigma_err,
    }
    return stats


def make_overlay_hist(values, hist_name, nbins, xmin, xmax):
    hist = ROOT.TH1D(hist_name, "", nbins, xmin, xmax)
    hist.Sumw2()
    for v in values:
        if math.isfinite(v):
            hist.Fill(v)
    return hist


def draw_distribution_overlay(
    sim_values, ref_values,
    sim_label, ref_label,
    title, x_title,
    out_png, hist_base_name,
    hist_range, nbins,
    normalize=True,
    text_pos="upper_left"
):
    xmin, xmax = hist_range

    h_sim = make_overlay_hist(sim_values, f"{hist_base_name}_sim", nbins, xmin, xmax)
    h_ref = make_overlay_hist(ref_values, f"{hist_base_name}_ref", nbins, xmin, xmax)

    sim_rows = len(sim_values)
    ref_rows = len(ref_values)

    if normalize:
        if h_sim.Integral() > 0:
            h_sim.Scale(1.0 / h_sim.Integral())
        if h_ref.Integral() > 0:
            h_ref.Scale(1.0 / h_ref.Integral())

    # REF: solid line
    h_ref.SetLineColor(ROOT.kBlue + 1)
    h_ref.SetLineWidth(3)
    h_ref.SetLineStyle(1)
    h_ref.SetFillStyle(0)

    # SIM: dashed line
    h_sim.SetLineColor(ROOT.kRed + 1)
    h_sim.SetLineWidth(3)
    h_sim.SetLineStyle(2)
    h_sim.SetFillStyle(0)

    c = ROOT.TCanvas(f"c_{hist_base_name}", "", 1000, 700)
    c.SetMargin(0.12, 0.04, 0.12, 0.06)

    ymax = max(h_sim.GetMaximum(), h_ref.GetMaximum()) * 1.35
    if ymax <= 0:
        ymax = 1.0

    h_ref.SetTitle(title)
    h_ref.GetXaxis().SetTitle(x_title)
    h_ref.GetYaxis().SetTitle("Normalized entries" if normalize else "Entries")
    h_ref.GetXaxis().CenterTitle()
    h_ref.GetYaxis().CenterTitle()
    h_ref.GetXaxis().SetTitleSize(0.045)
    h_ref.GetYaxis().SetTitleSize(0.045)
    h_ref.GetXaxis().SetLabelSize(0.04)
    h_ref.GetYaxis().SetLabelSize(0.04)
    h_ref.GetYaxis().SetMaxDigits(3)
    h_ref.SetMaximum(ymax)

    # draw reference first, then sim dashed on top
    h_ref.Draw("HIST")
    h_sim.Draw("HIST SAME")

    leg = ROOT.TLegend(0.72, 0.78, 0.97, 0.93)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(h_ref, ref_label, "l")
    leg.AddEntry(h_sim, sim_label, "l")
    leg.Draw()

    x0, y0, dy = choose_text_pos(text_pos)
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(dy)
    latex.DrawLatex(x0, y0,        f"{ref_label}: {ref_rows} rows")
    latex.DrawLatex(x0, y0 - 0.06, f"{sim_label}: {sim_rows} rows")

    c.SaveAs(out_png)

    return {"sim_hist": h_sim, "ref_hist": h_ref}


def compare_two(df_sim, df_ref, sim_label, ref_label, outdir, args, root_file):
    merge_keys = ["evt"]
    if args.merge_on == "evt_trk" and ("trkId" in df_sim.columns and "trkId" in df_ref.columns):
        merge_keys = ["evt", "trkId"]

    merged = pd.merge(
        df_sim[merge_keys + ["phi0", "omega"]].rename(columns={"phi0": "phi0_sim", "omega": "omega_sim"}),
        df_ref[merge_keys + ["phi0", "omega"]].rename(columns={"phi0": "phi0_ref", "omega": "omega_ref"}),
        on=merge_keys,
        how="inner"
    )

    if len(merged) == 0:
        raise RuntimeError("No matched rows between sim and ref")

    merged["dphi0_deg"] = merged.apply(
        lambda r: rad_to_deg(wrap_to_pi(float(r["phi0_sim"]) - float(r["phi0_ref"]))),
        axis=1
    )
    merged["domega"] = merged.apply(
        lambda r: float(r["omega_sim"]) - float(r["omega_ref"]),
        axis=1
    )

    merged["phi0_sim_deg"] = merged["phi0_sim"].apply(rad_to_deg)
    merged["phi0_ref_deg"] = merged["phi0_ref"].apply(rad_to_deg)

    safe = f"{sim_label}_minus_{ref_label}".replace(" ", "_")
    merged.to_csv(os.path.join(outdir, f"matched_{safe}.csv"), index=False)

    phi_stats = fit_and_draw(
        values=merged["dphi0_deg"].tolist(),
        title=f"phi0 Resolution: {sim_label} - {ref_label}",
        x_title=f"#Delta#phi_{{0}} = #phi_{{0}}^{{{sim_label}}} - #phi_{{0}}^{{{ref_label}}} [deg]",
        out_png=os.path.join(outdir, f"phi0_resolution_{safe}.png"),
        hist_name=f"hist_phi0_{safe}",
        hist_range=(args.phi0_xmin, args.phi0_xmax),
        fit_range=(args.phi0_fit_min, args.phi0_fit_max),
        nbins=args.phi0_bins,
        text_pos=args.phi0_text_pos,
    )

    omega_stats = fit_and_draw(
        values=merged["domega"].tolist(),
        title=f"omega Resolution: {sim_label} - {ref_label}",
        x_title=f"#Delta#omega = #omega^{{{sim_label}}} - #omega^{{{ref_label}}}",
        out_png=os.path.join(outdir, f"omega_resolution_{safe}.png"),
        hist_name=f"hist_omega_{safe}",
        hist_range=(args.omega_xmin, args.omega_xmax),
        fit_range=(args.omega_fit_min, args.omega_fit_max),
        nbins=args.omega_bins,
        text_pos=args.omega_text_pos,
    )

    phi_dist = draw_distribution_overlay(
        sim_values=merged["phi0_sim_deg"].tolist(),
        ref_values=merged["phi0_ref_deg"].tolist(),
        sim_label=sim_label,
        ref_label=ref_label,
        title=f"phi_{{0}} Distribution: {sim_label} vs {ref_label}",
        x_title="#phi_{0} [deg]",
        out_png=os.path.join(outdir, f"phi0_distribution_{safe}.png"),
        hist_base_name=f"dist_phi0_{safe}",
        hist_range=(args.phi0_dist_xmin, args.phi0_dist_xmax),
        nbins=args.phi0_dist_bins,
        normalize=args.normalize_distribution,
        text_pos=args.phi0_text_pos,
    )

    omega_dist = draw_distribution_overlay(
        sim_values=merged["omega_sim"].tolist(),
        ref_values=merged["omega_ref"].tolist(),
        sim_label=sim_label,
        ref_label=ref_label,
        title=f"omega Distribution: {sim_label} vs {ref_label}",
        x_title="#omega",
        out_png=os.path.join(outdir, f"omega_distribution_{safe}.png"),
        hist_base_name=f"dist_omega_{safe}",
        hist_range=(args.omega_dist_xmin, args.omega_dist_xmax),
        nbins=args.omega_dist_bins,
        normalize=args.normalize_distribution,
        text_pos=args.omega_text_pos,
    )

    root_file.cd()
    phi_stats["hist"].Write(f"phi0_diff_hist_{safe}")
    phi_stats["graph"].Write(f"phi0_diff_graph_{safe}")
    phi_stats["fit"].Write(f"phi0_diff_fit_{safe}")

    omega_stats["hist"].Write(f"omega_diff_hist_{safe}")
    omega_stats["graph"].Write(f"omega_diff_graph_{safe}")
    omega_stats["fit"].Write(f"omega_diff_fit_{safe}")

    phi_dist["sim_hist"].Write(f"phi0_dist_sim_{safe}")
    phi_dist["ref_hist"].Write(f"phi0_dist_ref_{safe}")
    omega_dist["sim_hist"].Write(f"omega_dist_sim_{safe}")
    omega_dist["ref_hist"].Write(f"omega_dist_ref_{safe}")

    return merged, phi_stats, omega_stats


def print_stats(label, stats, unit_suffix=""):
    print(f"[{label}]")
    print(f"  Entries          = {stats['entries']:.0f}")
    print(f"  Hist mean        = {stats['hist_mean']:.10g} +/- {stats['hist_mean_err']:.4g}{unit_suffix}")
    print(f"  Hist std         = {stats['hist_std']:.10g} +/- {stats['hist_std_err']:.4g}{unit_suffix}")
    if math.isfinite(stats["fit_mean"]):
        print(f"  Fit mean         = {stats['fit_mean']:.10g} +/- {stats['fit_mean_err']:.4g}{unit_suffix}")
        print(f"  Fit sigma        = {stats['fit_sigma']:.10g} +/- {stats['fit_sigma_err']:.4g}{unit_suffix}")
    else:
        print("  Fit failed or not available")


def main():
    parser = argparse.ArgumentParser(
        description="Compare simulation CSV against reference CSV with resolution and distribution plots"
    )
    parser.add_argument("--sim", required=True, help="CSV for simulation data")
    parser.add_argument("--ref", required=True, help="CSV for reference data")
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--merge-on", choices=["evt", "evt_trk"], default="evt",
                        help="Use evt or evt+trkId for matching")
    parser.add_argument("--sim-trkid", type=int, default=None,
                        help="Optional: keep only rows with this trkId in sim CSV before comparison")
    parser.add_argument("--ref-trkid", type=int, default=None,
                        help="Optional: keep only rows with this trkId in ref CSV before comparison")

    parser.add_argument("--label-sim", default="SIM")
    parser.add_argument("--label-ref", default="REF")

    # resolution ranges
    parser.add_argument("--phi0-xmin", type=float, default=-2.0)
    parser.add_argument("--phi0-xmax", type=float, default= 2.0)
    parser.add_argument("--phi0-fit-min", type=float, default=-2.0)
    parser.add_argument("--phi0-fit-max", type=float, default= 2.0)
    parser.add_argument("--phi0-bins", type=int, default=40)

    parser.add_argument("--omega-xmin", type=float, default=-0.0015)
    parser.add_argument("--omega-xmax", type=float, default= 0.0015)
    parser.add_argument("--omega-fit-min", type=float, default=-0.0010)
    parser.add_argument("--omega-fit-max", type=float, default= 0.0010)
    parser.add_argument("--omega-bins", type=int, default=50)

    # distribution ranges
    parser.add_argument("--phi0-dist-xmin", type=float, default=-180.0)
    parser.add_argument("--phi0-dist-xmax", type=float, default= 180.0)
    parser.add_argument("--phi0-dist-bins", type=int, default=72)

    parser.add_argument("--omega-dist-xmin", type=float, default=-0.03)
    parser.add_argument("--omega-dist-xmax", type=float, default= 0.03)
    parser.add_argument("--omega-dist-bins", type=int, default=80)

    parser.add_argument("--normalize-distribution", action="store_true",
                        help="Normalize distribution overlays")

    parser.add_argument("--phi0-text-pos", choices=["upper_right_far", "upper_right", "upper_left"],
                        default="upper_left")
    parser.add_argument("--omega-text-pos", choices=["upper_right_far", "upper_right", "upper_left"],
                        default="upper_left")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    sim = pd.read_csv(args.sim)
    ref = pd.read_csv(args.ref)

    required_cols = {"evt", "phi0", "omega"}
    for name, df in [("SIM", sim), ("REF", ref)]:
        if not required_cols.issubset(df.columns):
            raise RuntimeError(f"{name} CSV must contain evt, phi0, omega columns")

    if args.sim_trkid is not None:
        if "trkId" not in sim.columns:
            raise RuntimeError("--sim-trkid was set but sim CSV has no trkId column")
        sim = sim[sim["trkId"] == args.sim_trkid].copy()
        sim.to_csv(os.path.join(args.outdir, "filtered_sim.csv"), index=False)

    if args.ref_trkid is not None:
        if "trkId" not in ref.columns:
            raise RuntimeError("--ref-trkid was set but ref CSV has no trkId column")
        ref = ref[ref["trkId"] == args.ref_trkid].copy()
        ref.to_csv(os.path.join(args.outdir, "filtered_ref.csv"), index=False)

    root_out = ROOT.TFile.Open(os.path.join(args.outdir, "comparison.root"), "RECREATE")

    merged, phi_stats, omega_stats = compare_two(
        sim, ref,
        args.label_sim, args.label_ref,
        args.outdir, args, root_out
    )

    root_out.Write()
    root_out.Close()

    print(f"SIM rows used: {len(sim)}")
    print(f"REF rows used: {len(ref)}")
    print(f"Matched rows: {len(merged)}")
    print_stats("phi0 Resolution", phi_stats, " deg")
    print_stats("omega Resolution", omega_stats, "")
    print(f"\nOutputs written to: {args.outdir}")


if __name__ == "__main__":
    main()