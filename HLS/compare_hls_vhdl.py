#!/usr/bin/env python3
"""
compare_hls_vhdl.py
====================
Compare two CSV files (e.g. HLS C-sim vs Vivado VHDL sim) on rho and phi0.

Expected CSV columns:
    evt, trkId, phi0 [rad], rho [cm]
    (T0 and idx columns are ignored if present)

Usage examples
--------------
# HLS vs VHDL:
python compare_hls_vhdl.py \\
    --sim hls_local_aligned_out_nts.csv \\
    --ref vivado_out.csv \\
    --label-sim HLS --label-ref VHDL \\
    --outdir ./results_hls_vs_vhdl

# HLS vs MCP (if you have a MCP reference CSV with evt,trkId,phi0,rho):
python compare_hls_vhdl.py \\
    --sim hls_local_aligned_out_nts.csv \\
    --ref mcp_reference.csv \\
    --label-sim HLS --label-ref MCP \\
    --outdir ./results_hls_vs_mcp

Outputs (all in --outdir):
    matched_<SIM>_minus_<REF>.csv         merged rows with diff columns
    phi0_resolution_<SIM>_minus_<REF>.png  Δphi0 histogram + Gaussian fit
    rho_resolution_<SIM>_minus_<REF>.png   Δrho  histogram + Gaussian fit
    phi0_distribution_<SIM>_vs_<REF>.png   phi0 overlay (normalized)
    rho_distribution_<SIM>_vs_<REF>.png    rho  overlay (normalized)
    comparison.root                         all histograms + fits
    summary.txt                             printed stats
"""

import os
import math
import argparse
import pandas as pd
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)


# ── Utilities ──────────────────────────────────────────────────────────────────

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
    return 0.16, 0.88, 0.040   # upper_left (default)


# ── Single-variable resolution plot ───────────────────────────────────────────

def make_graph_from_hist(hist):
    """Build TGraphAsymmErrors with Poisson error bars from a TH1D."""
    xv, yv, exl, exh, eyl, eyh = [], [], [], [], [], []
    for ib in range(1, hist.GetNbinsX() + 1):
        y = hist.GetBinContent(ib)
        if y <= 0:
            continue
        x   = hist.GetBinCenter(ib)
        bw  = 0.5 * hist.GetBinWidth(ib)
        err = math.sqrt(y)
        xv.append(x);  yv.append(y)
        exl.append(bw); exh.append(bw)
        eyl.append(err); eyh.append(err)
    g = ROOT.TGraphAsymmErrors(
        len(xv),
        array("d", xv), array("d", yv),
        array("d", exl), array("d", exh),
        array("d", eyl), array("d", eyh),
    )
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    g.SetLineWidth(2)
    return g


def resolution_plot(values, title, x_title, out_png, hist_name,
                    hist_range, fit_range, nbins, text_pos,
                    fit_model="gaus"):
    """
    Draw a resolution histogram with optional Gaussian or Double-Gaussian fit.
    fit_model: "gaus" | "dgaus" | "none"
    Returns stats dict.
    """
    xmin, xmax = hist_range
    fmin, fmax = fit_range

    hist = ROOT.TH1D(hist_name, "", nbins, xmin, xmax)
    hist.Sumw2()
    for v in values:
        if math.isfinite(v):
            hist.Fill(v)

    graph = make_graph_from_hist(hist)

    c = ROOT.TCanvas(f"c_{hist_name}", "", 900, 700)
    c.SetMargin(0.12, 0.04, 0.12, 0.06)

    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle("Entries")
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.040)
    hist.GetYaxis().SetLabelSize(0.040)
    hist.GetYaxis().SetMaxDigits(3)
    hist.SetMaximum(max(hist.GetMaximum() * 1.45, 1.0))
    hist.Draw("HIST")
    graph.Draw("P E SAME")

    mean_h    = hist.GetMean()
    mean_h_e  = hist.GetMeanError()
    std_h     = hist.GetStdDev()
    std_h_e   = hist.GetStdDevError()

    fit_results = {}

    if fit_model == "gaus":
        fit = ROOT.TF1(f"fit_{hist_name}", "gaus", fmin, fmax)
        fit.SetLineWidth(2)
        fit.SetLineColor(ROOT.kRed + 1)
        r = hist.Fit(fit, "RQ0S")
        if r and int(r) == 0:
            fit.Draw("SAME")
            fit_results = {
                "mean":  fit.GetParameter(1), "mean_e":  fit.GetParError(1),
                "sigma": fit.GetParameter(2), "sigma_e": fit.GetParError(2),
            }

    elif fit_model == "dgaus":
        # Double Gaussian:  N1*gaus(m,s1) + N2*gaus(m,s2)
        dgaus_str = ("[0]*exp(-0.5*((x-[1])/[2])^2)"
                     "+[3]*exp(-0.5*((x-[1])/[4])^2)")
        fit = ROOT.TF1(f"fit_{hist_name}", dgaus_str, fmin, fmax)
        amp = hist.GetMaximum()
        fit.SetParameters(amp * 0.4, 0.0, std_h * 0.3,
                          amp * 0.6, std_h)
        fit.SetParLimits(2, 1e-9, abs(fmax - fmin))
        fit.SetParLimits(4, 1e-9, abs(fmax - fmin))
        fit.SetLineWidth(2)
        fit.SetLineColor(ROOT.kRed + 1)
        r = hist.Fit(fit, "RQ0S")
        if r and int(r) == 0:
            fit.Draw("SAME")
            s1  = abs(fit.GetParameter(2))
            s2  = abs(fit.GetParameter(4))
            A1  = abs(fit.GetParameter(0))
            A2  = abs(fit.GetParameter(3))
            f1  = A1 / (A1 + A2)
            f2  = 1.0 - f1
            if s1 <= s2:
                core_s, tail_s, core_f, tail_f = s1, s2, f1, f2
            else:
                core_s, tail_s, core_f, tail_f = s2, s1, f2, f1
            sig_eff = math.sqrt(f1 * s1**2 + f2 * s2**2)
            fit_results = {
                "mean":       fit.GetParameter(1),
                "mean_e":     fit.GetParError(1),
                "core_sigma": core_s,
                "tail_sigma": tail_s,
                "core_frac":  core_f,
                "tail_frac":  tail_f,
                "sigma_eff":  sig_eff,
            }

    # ── Text overlay ──────────────────────────────────────────────────────
    x0, y0, dy = choose_text_pos(text_pos)
    lat = ROOT.TLatex()
    lat.SetNDC(True)
    lat.SetTextSize(dy)

    n_valid = sum(1 for v in values if math.isfinite(v))
    n_total = len(values)

    lat.DrawLatex(x0, y0,          f"Total entries   = {n_total}")
    lat.DrawLatex(x0, y0 - 0.06,   f"Valid entries   = {n_valid}")
    lat.DrawLatex(x0, y0 - 0.12,   f"Hist mean = {mean_h:.6g} #pm {mean_h_e:.3g}")
    lat.DrawLatex(x0, y0 - 0.18,   f"Hist std  = {std_h:.6g} #pm {std_h_e:.3g}")

    row = 0.24
    if fit_model == "gaus" and fit_results:
        lat.DrawLatex(x0, y0 - row,
            f"Fit mean  = {fit_results['mean']:.6g} #pm {fit_results['mean_e']:.3g}")
        lat.DrawLatex(x0, y0 - row - 0.06,
            f"Fit #sigma = {fit_results['sigma']:.6g} #pm {fit_results['sigma_e']:.3g}")

    elif fit_model == "dgaus" and fit_results:
        lat.DrawLatex(x0, y0 - row,
            f"Fit mean     = {fit_results['mean']:.6g} #pm {fit_results['mean_e']:.3g}")
        lat.DrawLatex(x0, y0 - row - 0.06,
            f"Core #sigma  = {fit_results['core_sigma']:.6g}")
        lat.DrawLatex(x0, y0 - row - 0.12,
            f"Tail #sigma  = {fit_results['tail_sigma']:.6g}")
        lat.DrawLatex(x0, y0 - row - 0.18,
            f"#sigma_{{eff}} = {fit_results['sigma_eff']:.6g}")

    c.SaveAs(out_png)
    c.Close()

    return {
        "hist": hist, "graph": graph,
        "n_total": n_total, "n_valid": n_valid,
        "hist_mean": mean_h, "hist_mean_e": mean_h_e,
        "hist_std":  std_h,  "hist_std_e":  std_h_e,
        **fit_results,
    }


# ── Distribution overlay ───────────────────────────────────────────────────────

def distribution_overlay(sim_vals, ref_vals,
                          sim_label, ref_label,
                          title, x_title,
                          out_png, hist_base,
                          hist_range, nbins,
                          normalize=True, text_pos="upper_left"):
    xmin, xmax = hist_range

    def fill(name, vals):
        h = ROOT.TH1D(name, "", nbins, xmin, xmax)
        h.Sumw2()
        for v in vals:
            if math.isfinite(v):
                h.Fill(v)
        return h

    h_ref = fill(f"{hist_base}_ref", ref_vals)
    h_sim = fill(f"{hist_base}_sim", sim_vals)

    if normalize:
        for h in (h_ref, h_sim):
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

    h_ref.SetLineColor(ROOT.kBlue + 1); h_ref.SetLineWidth(3); h_ref.SetLineStyle(1); h_ref.SetFillStyle(0)
    h_sim.SetLineColor(ROOT.kRed  + 1); h_sim.SetLineWidth(3); h_sim.SetLineStyle(2); h_sim.SetFillStyle(0)

    c = ROOT.TCanvas(f"c_{hist_base}", "", 1000, 700)
    c.SetMargin(0.12, 0.04, 0.12, 0.06)

    ymax = max(h_ref.GetMaximum(), h_sim.GetMaximum()) * 1.40
    if ymax <= 0:
        ymax = 1.0

    h_ref.SetTitle(title)
    h_ref.GetXaxis().SetTitle(x_title)
    h_ref.GetYaxis().SetTitle("Normalized entries" if normalize else "Entries")
    h_ref.GetXaxis().CenterTitle(); h_ref.GetYaxis().CenterTitle()
    h_ref.GetXaxis().SetTitleSize(0.045); h_ref.GetYaxis().SetTitleSize(0.045)
    h_ref.GetXaxis().SetLabelSize(0.040); h_ref.GetYaxis().SetLabelSize(0.040)
    h_ref.GetYaxis().SetMaxDigits(3)
    h_ref.SetMaximum(ymax)
    h_ref.Draw("HIST")
    h_sim.Draw("HIST SAME")

    leg = ROOT.TLegend(0.72, 0.80, 0.97, 0.93)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_ref, ref_label, "l")
    leg.AddEntry(h_sim, sim_label, "l")
    leg.Draw()

    x0, y0, dy = choose_text_pos(text_pos)
    lat = ROOT.TLatex(); lat.SetNDC(True); lat.SetTextSize(dy)
    lat.DrawLatex(x0, y0,        f"{ref_label}: {len(ref_vals)} tracks")
    lat.DrawLatex(x0, y0 - 0.06, f"{sim_label}: {len(sim_vals)} tracks")

    c.SaveAs(out_png)
    c.Close()

    return {"sim_hist": h_sim, "ref_hist": h_ref}


# ── Main comparison ────────────────────────────────────────────────────────────

def compare(df_sim, df_ref, sim_label, ref_label, outdir, args, root_file):
    # Merge on evt (+ trkId if available in both)
    merge_keys = ["evt"]
    if ("trkId" in df_sim.columns and "trkId" in df_ref.columns
            and args.merge_on == "evt_trk"):
        merge_keys = ["evt", "trkId"]

    merged = pd.merge(
        df_sim[merge_keys + ["phi0", "rho"]].rename(
            columns={"phi0": "phi0_sim", "rho": "rho_sim"}),
        df_ref[merge_keys + ["phi0", "rho"]].rename(
            columns={"phi0": "phi0_ref", "rho": "rho_ref"}),
        on=merge_keys, how="inner",
    )

    if len(merged) == 0:
        raise RuntimeError("No matched rows after merge — check evt column values")

    safe = f"{sim_label}_minus_{ref_label}".replace(" ", "_")

    # ── Differences ──────────────────────────────────────────────────────
    merged["dphi0_deg"] = merged.apply(
        lambda r: rad_to_deg(wrap_to_pi(float(r["phi0_sim"]) - float(r["phi0_ref"]))),
        axis=1,
    )
    merged["drho"] = merged["rho_sim"].astype(float) - merged["rho_ref"].astype(float)
    merged["drho_rel"] = merged.apply(
        lambda r: (float(r["rho_sim"]) - float(r["rho_ref"])) / float(r["rho_ref"])
                  if float(r["rho_ref"]) != 0 else float("nan"),
        axis=1,
    )
    merged["phi0_sim_deg"] = merged["phi0_sim"].apply(lambda x: rad_to_deg(float(x)))
    merged["phi0_ref_deg"] = merged["phi0_ref"].apply(lambda x: rad_to_deg(float(x)))

    merged.to_csv(os.path.join(outdir, f"matched_{safe}.csv"), index=False)

    # ── phi0 resolution ──────────────────────────────────────────────────
    phi0_stats = resolution_plot(
        values    = merged["dphi0_deg"].tolist(),
        title     = f"#phi_{{0}} Resolution: {sim_label} #minus {ref_label}",
        x_title   = (f"#Delta#phi_{{0}} = #phi_{{0}}^{{{sim_label}}}"
                     f" #minus #phi_{{0}}^{{{ref_label}}} [deg]"),
        out_png   = os.path.join(outdir, f"phi0_resolution_{safe}.png"),
        hist_name = f"h_phi0_res_{safe}",
        hist_range= (args.phi0_xmin, args.phi0_xmax),
        fit_range = (args.phi0_fit_min, args.phi0_fit_max),
        nbins     = args.phi0_bins,
        text_pos  = args.phi0_text_pos,
        fit_model = args.phi0_fit,
    )

    # ── rho resolution (absolute) ────────────────────────────────────────
    rho_stats = resolution_plot(
        values    = merged["drho"].tolist(),
        title     = f"#rho Resolution: {sim_label} #minus {ref_label}",
        x_title   = (f"#Delta#rho = #rho^{{{sim_label}}}"
                     f" #minus #rho^{{{ref_label}}} [cm]"),
        out_png   = os.path.join(outdir, f"rho_resolution_{safe}.png"),
        hist_name = f"h_rho_res_{safe}",
        hist_range= (args.rho_xmin, args.rho_xmax),
        fit_range = (args.rho_fit_min, args.rho_fit_max),
        nbins     = args.rho_bins,
        text_pos  = args.rho_text_pos,
        fit_model = args.rho_fit,
    )

    # ── rho resolution (relative %) ──────────────────────────────────────
    rho_rel_stats = resolution_plot(
        values    = [v * 100.0 for v in merged["drho_rel"].tolist()],
        title     = f"#rho Relative Resolution: {sim_label} #minus {ref_label}",
        x_title   = (f"#Delta#rho/#rho [%]"),
        out_png   = os.path.join(outdir, f"rho_rel_resolution_{safe}.png"),
        hist_name = f"h_rho_rel_res_{safe}",
        hist_range= (args.rho_rel_xmin, args.rho_rel_xmax),
        fit_range = (args.rho_rel_fit_min, args.rho_rel_fit_max),
        nbins     = args.rho_rel_bins,
        text_pos  = args.rho_text_pos,
        fit_model = args.rho_fit,
    )

    # ── phi0 distribution overlay ────────────────────────────────────────
    phi0_dist = distribution_overlay(
        sim_vals  = merged["phi0_sim_deg"].tolist(),
        ref_vals  = merged["phi0_ref_deg"].tolist(),
        sim_label = sim_label, ref_label = ref_label,
        title     = f"#phi_{{0}} Distribution: {sim_label} vs {ref_label}",
        x_title   = "#phi_{0} [deg]",
        out_png   = os.path.join(outdir, f"phi0_distribution_{safe}.png"),
        hist_base = f"dist_phi0_{safe}",
        hist_range= (args.phi0_dist_xmin, args.phi0_dist_xmax),
        nbins     = args.phi0_dist_bins,
        normalize = args.normalize_distribution,
        text_pos  = args.phi0_text_pos,
    )

    # ── rho distribution overlay ─────────────────────────────────────────
    rho_dist = distribution_overlay(
        sim_vals  = merged["rho_sim"].astype(float).tolist(),
        ref_vals  = merged["rho_ref"].astype(float).tolist(),
        sim_label = sim_label, ref_label = ref_label,
        title     = f"#rho Distribution: {sim_label} vs {ref_label}",
        x_title   = "#rho [cm]",
        out_png   = os.path.join(outdir, f"rho_distribution_{safe}.png"),
        hist_base = f"dist_rho_{safe}",
        hist_range= (args.rho_dist_xmin, args.rho_dist_xmax),
        nbins     = args.rho_dist_bins,
        normalize = args.normalize_distribution,
        text_pos  = args.rho_text_pos,
    )

    # ── Save to ROOT file ─────────────────────────────────────────────────
    root_file.cd()
    for name, obj in [
        (f"phi0_res_hist_{safe}",   phi0_stats["hist"]),
        (f"phi0_res_graph_{safe}",  phi0_stats["graph"]),
        (f"rho_res_hist_{safe}",    rho_stats["hist"]),
        (f"rho_res_graph_{safe}",   rho_stats["graph"]),
        (f"rho_rel_hist_{safe}",    rho_rel_stats["hist"]),
        (f"phi0_dist_sim_{safe}",   phi0_dist["sim_hist"]),
        (f"phi0_dist_ref_{safe}",   phi0_dist["ref_hist"]),
        (f"rho_dist_sim_{safe}",    rho_dist["sim_hist"]),
        (f"rho_dist_ref_{safe}",    rho_dist["ref_hist"]),
    ]:
        obj.Write(name)

    return merged, phi0_stats, rho_stats, rho_rel_stats


# ── Print summary ──────────────────────────────────────────────────────────────

def print_stats(label, stats, unit=""):
    lines = [
        f"[{label}]",
        f"  Total entries   = {stats['n_total']}",
        f"  Valid entries   = {stats['n_valid']}",
        f"  Hist mean       = {stats['hist_mean']:.10g} +/- {stats['hist_mean_e']:.4g}{unit}",
        f"  Hist std        = {stats['hist_std']:.10g} +/- {stats['hist_std_e']:.4g}{unit}",
    ]
    # Gaussian fit
    if "sigma" in stats and math.isfinite(stats.get("sigma", float("nan"))):
        lines += [
            f"  Fit mean        = {stats['mean']:.10g} +/- {stats['mean_e']:.4g}{unit}",
            f"  Fit sigma       = {stats['sigma']:.10g} +/- {stats['sigma_e']:.4g}{unit}",
        ]
    # Double-Gaussian fit
    if "core_sigma" in stats:
        lines += [
            f"  Fit mean        = {stats['mean']:.10g} +/- {stats['mean_e']:.4g}{unit}",
            f"  Core sigma      = {stats['core_sigma']:.10g}{unit}",
            f"  Tail sigma      = {stats['tail_sigma']:.10g}{unit}",
            f"  Core fraction   = {stats['core_frac']:.4f}",
            f"  sigma_eff       = {stats['sigma_eff']:.10g}{unit}",
        ]
    return "\n".join(lines)


# ── CLI ───────────────────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--sim",       required=True,  help="Simulation CSV (HLS csim or RTL sim)")
    p.add_argument("--ref",       required=True,  help="Reference CSV (VHDL or MCP)")
    p.add_argument("--outdir",    required=True,  help="Output directory")
    p.add_argument("--label-sim", default="HLS",  help="Label for sim data")
    p.add_argument("--label-ref", default="VHDL", help="Label for ref data")
    p.add_argument("--merge-on",  choices=["evt","evt_trk"], default="evt_trk",
                   help="Merge key (default: evt+trkId)")

    # Filter by trkId before merge
    p.add_argument("--sim-trkid", type=int, default=None,
                   help="Keep only this trkId in sim CSV before comparison")
    p.add_argument("--ref-trkid", type=int, default=None,
                   help="Keep only this trkId in ref CSV before comparison")

    # phi0 resolution
    p.add_argument("--phi0-xmin",     type=float, default=-2.0)
    p.add_argument("--phi0-xmax",     type=float, default= 2.0)
    p.add_argument("--phi0-fit-min",  type=float, default=-2.0)
    p.add_argument("--phi0-fit-max",  type=float, default= 2.0)
    p.add_argument("--phi0-bins",     type=int,   default=80)
    p.add_argument("--phi0-fit",      choices=["gaus","dgaus","none"], default="dgaus",
                   help="Fit model for phi0 resolution")
    p.add_argument("--phi0-text-pos", choices=["upper_right_far","upper_right","upper_left"],
                   default="upper_left")

    # rho resolution (absolute)
    p.add_argument("--rho-xmin",     type=float, default=-5.0,
                   help="Δrho axis min [cm]")
    p.add_argument("--rho-xmax",     type=float, default= 5.0)
    p.add_argument("--rho-fit-min",  type=float, default=-5.0)
    p.add_argument("--rho-fit-max",  type=float, default= 5.0)
    p.add_argument("--rho-bins",     type=int,   default=80)
    p.add_argument("--rho-fit",      choices=["gaus","dgaus","none"], default="gaus",
                   help="Fit model for rho resolution")
    p.add_argument("--rho-text-pos", choices=["upper_right_far","upper_right","upper_left"],
                   default="upper_left")

    # rho relative resolution
    p.add_argument("--rho-rel-xmin",    type=float, default=-1.0,
                   help="Δrho/rho axis min [%%]")
    p.add_argument("--rho-rel-xmax",    type=float, default= 1.0)
    p.add_argument("--rho-rel-fit-min", type=float, default=-1.0)
    p.add_argument("--rho-rel-fit-max", type=float, default= 1.0)
    p.add_argument("--rho-rel-bins",    type=int,   default=80)

    # phi0 distribution
    p.add_argument("--phi0-dist-xmin", type=float, default=-180.0)
    p.add_argument("--phi0-dist-xmax", type=float, default= 180.0)
    p.add_argument("--phi0-dist-bins", type=int,   default=72)

    # rho distribution
    p.add_argument("--rho-dist-xmin",  type=float, default=0.0)
    p.add_argument("--rho-dist-xmax",  type=float, default=2000.0)
    p.add_argument("--rho-dist-bins",  type=int,   default=80)

    p.add_argument("--normalize-distribution", action="store_true",
                   help="Normalize distribution overlays to unit area")

    return p


def main():
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load CSVs ─────────────────────────────────────────────────────────
    sim = pd.read_csv(args.sim)
    ref = pd.read_csv(args.ref)

    for name, df in [("SIM", sim), ("REF", ref)]:
        for col in ("evt", "phi0", "rho"):
            if col not in df.columns:
                raise RuntimeError(f"{name} CSV is missing column '{col}'")

    if args.sim_trkid is not None:
        sim = sim[sim["trkId"] == args.sim_trkid].copy()
    if args.ref_trkid is not None:
        ref = ref[ref["trkId"] == args.ref_trkid].copy()

    print(f"SIM rows after filter: {len(sim)}")
    print(f"REF rows after filter: {len(ref)}")

    root_out = ROOT.TFile.Open(
        os.path.join(args.outdir, "comparison.root"), "RECREATE"
    )

    merged, phi0_stats, rho_stats, rho_rel_stats = compare(
        sim, ref,
        args.label_sim, args.label_ref,
        args.outdir, args, root_out,
    )

    root_out.Write()
    root_out.Close()

    # ── Print summary ─────────────────────────────────────────────────────
    summary_lines = [
        f"SIM: {args.sim}",
        f"REF: {args.ref}",
        f"Matched rows: {len(merged)}",
        "",
        print_stats(f"phi0 Resolution [{args.label_sim} - {args.label_ref}]",
                    phi0_stats, " deg"),
        "",
        print_stats(f"rho Resolution [{args.label_sim} - {args.label_ref}]",
                    rho_stats, " cm"),
        "",
        print_stats(f"rho Relative Resolution [{args.label_sim} - {args.label_ref}]",
                    rho_rel_stats, " %"),
        "",
        f"Outputs: {args.outdir}",
    ]

    summary = "\n".join(summary_lines)
    print(summary)

    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(summary + "\n")


if __name__ == "__main__":
    main()