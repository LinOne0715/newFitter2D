#include "cdc2dfit_kernel.h"
#include "lut.h"
#include "drift_phi_lut_generated.h"
#include "sincos_lut_generated.h"
#include "atan_unit_lut_generated.h"
#include "sqrt_lut_generated.h"

static inline ap_int<4> sl_raw_to_ax(ap_uint<4> sl_raw) {
#pragma HLS INLINE
    return ap_int<4>(sl_raw >> 1);
}

static inline fixed_type_general my_atan2_with_recip(
    fixed_type_general y,  fixed_type_general x,
    fixed_type_general recip_ax, fixed_type_general recip_ay
) {
#pragma HLS INLINE
    fixed_type_general ax = (x < FX_ZERO) ? fixed_type_general(-x) : x;
    fixed_type_general ay = (y < FX_ZERO) ? fixed_type_general(-y) : y;
    fixed_type_general base;
    if (ax >= ay) {
        fixed_type_general z = (ax != FX_ZERO) ? fixed_type_general(ay * recip_ax) : FX_ZERO;
        base = call_atan_unit_lut(z);
    } else {
        fixed_type_general z = (ay != FX_ZERO) ? fixed_type_general(ax * recip_ay) : FX_ZERO;
        base = FX_HALF_PI - call_atan_unit_lut(z);
    }
    if (x >= FX_ZERO) return (y >= FX_ZERO) ? base : fixed_type_general(-base);
    return (y >= FX_ZERO) ? fixed_type_general(FX_PI - base) : fixed_type_general(base - FX_PI);
}

static inline fixed_type_general wrap_to_npi_pi(fixed_type_general phi) {
#pragma HLS INLINE
    if      (phi >  FX_PI)  phi -= FX_TWO_PI;
    else if (phi <= -FX_PI) phi += FX_TWO_PI;
    return phi;
}

static inline fixed_type_general stored_track_phi0(
    fixed_type_general phi0_internal, ap_int<3> chargeFit) {
#pragma HLS INLINE
    // chargeFit in {-1,0,1}: 3-bit × 32-bit = 1×2 = 2 DSP (was 3×64 = 4 DSP)
    typedef ap_fixed<32,4> half_pi_t;
    return wrap_to_npi_pi(phi0_internal + fixed_type_general(
        fixed_type_general(chargeFit) * half_pi_t(FX_HALF_PI)));
}

static inline fixed_type_general exported_phi0_like_basf2(fixed_type_general phi0_internal) {
#pragma HLS INLINE
    return wrap_to_npi_pi(phi0_internal);
}

static inline fixed_type_general calPhi_basf2_lut(
    fixed_type_general wirePhi, ap_int<4> ax, int dt, ap_uint<2> lr) {
#pragma HLS INLINE
    fixed_type_general driftPhi = call_drift_phi_lut(ax, dt);
    if      (lr == 1) return fixed_type_general(wirePhi - driftPhi);
    else if (lr == 2) return fixed_type_general(wirePhi + driftPhi);
    return wirePhi;
}

static inline void process_hit_local_aligned(
    ap_uint<6> layer, ap_uint<4> sl_raw, ap_uint<10> wire, ap_uint<2> lr,
    ap_uint<10> priorityTime, ap_uint<14> driftTick, ap_int<16> T0,
    fixed_type_general &phi, fixed_type_general &invErr,
    fixed_type_general &tsId, fixed_type_general &useSL, ap_uint<1> &used
) {
#pragma HLS INLINE
    (void)priorityTime;
    ap_int<4> ax = sl_raw_to_ax(sl_raw);
    if (ax < 0 || ax > 4 || layer != PRIORITY_LAYERS[(int)ax]) {
        phi = invErr = tsId = useSL = FX_ZERO; used = 0; return;
    }
    int dt = (int)driftTick;
    if (dt < 0) dt = 0; if (dt > 511) dt = 511;
    fixed_type_general wirePhi = fixed_type_general(wire) * WIRE_STEP[(int)ax];
    phi    = (T0 == 9999) ? wirePhi : calPhi_basf2_lut(wirePhi, ax, dt, lr);
    invErr = ((lr != 3) && (T0 != 9999)) ? INV_DRIFT_PHI_ERR[(int)ax] : INV_WIRE_PHI_ERR[(int)ax];
    tsId = fixed_type_general(wire); useSL = FX_ONE; used = 1;
}

// ── Stage 1 ──────────────────────────────────────────────────────────────────
static void stage1_hit(
    ap_uint<6>  layer[5], ap_uint<4>  sl_raw[5], ap_uint<10> wire[5],
    ap_uint<2>  lr[5],    ap_uint<10> pTime[5],  ap_uint<14> dTick[5],
    ap_int<16>  T0,
    fixed_type_general phi_o[5],    fixed_type_general invErr_o[5],
    fixed_type_general tsId_o[5],   fixed_type_general useSL_o[5],
    ap_uint<3> &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=layer    complete
#pragma HLS ARRAY_PARTITION variable=sl_raw   complete
#pragma HLS ARRAY_PARTITION variable=wire     complete
#pragma HLS ARRAY_PARTITION variable=lr       complete
#pragma HLS ARRAY_PARTITION variable=pTime    complete
#pragma HLS ARRAY_PARTITION variable=dTick    complete
#pragma HLS ARRAY_PARTITION variable=phi_o    complete
#pragma HLS ARRAY_PARTITION variable=invErr_o complete
#pragma HLS ARRAY_PARTITION variable=tsId_o   complete
#pragma HLS ARRAY_PARTITION variable=useSL_o  complete

    ap_uint<1> used[5];
#pragma HLS ARRAY_PARTITION variable=used complete
    for (int i = 0; i < 5; ++i) {
#pragma HLS UNROLL
        process_hit_local_aligned(
            layer[i], sl_raw[i], wire[i], lr[i], pTime[i], dTick[i], T0,
            phi_o[i], invErr_o[i], tsId_o[i], useSL_o[i], used[i]);
    }
    nhits_o = used[0] + used[1] + used[2] + used[3] + used[4];
}

// ── Stage 2a: cos/sin + accumulate A..E ──────────────────────────────────────
// ROM_NP sincos + 5 UNROLL = parallel reads -> II=1
static void stage2a_accum(
    fixed_type_general phi[5],   fixed_type_general invErr[5],
    ap_uint<3>  nhits_i,
    fixed_type_accum_ABC &A_o, fixed_type_accum_ABC &B_o, fixed_type_accum_ABC &C_o,
    fixed_type_accum_DE  &D_o, fixed_type_accum_DE  &E_o,
    ap_uint<3>  &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=phi    complete
#pragma HLS ARRAY_PARTITION variable=invErr complete

    nhits_o = nhits_i;
    A_o = B_o = C_o = fixed_type_accum_ABC(0);
    D_o = E_o = fixed_type_accum_DE(0);
    if (nhits_i < 2) return;

    fixed_type_radius r[5] = {
        (fixed_type_radius)RADIUS[0], (fixed_type_radius)RADIUS[1],
        (fixed_type_radius)RADIUS[2], (fixed_type_radius)RADIUS[3],
        (fixed_type_radius)RADIUS[4]
    };
#pragma HLS ARRAY_PARTITION variable=r complete

    fixed_type_invF invF[5];
#pragma HLS ARRAY_PARTITION variable=invF complete

    for (int i = 0; i < 5; ++i) {
#pragma HLS UNROLL
        // invErr(24bit) * INV_RADIUS(14bit) = 24×14 -> 1 DSP
        fixed_type_const ie_n = (fixed_type_const)invErr[i];
        fixed_type_const ir_n = (fixed_type_const)INV_RADIUS[i];
        invF[i] = (invErr[i] == FX_ZERO) ? fixed_type_invF(0) : fixed_type_invF(ie_n * ir_n);
    }

    for (int i = 0; i < 5; ++i) {
#pragma HLS UNROLL
        // w: invF(17) * invF(17) = 17×17 -> 1 DSP
        fixed_type_weight w  = fixed_type_weight(invF[i] * invF[i]);
        // cos/sin: 20-bit trig
        fixed_type_trig   c  = (fixed_type_trig)call_cos_lut(phi[i]);
        fixed_type_trig   s  = (fixed_type_trig)call_sin_lut(phi[i]);
        // cw,sw: trig(20) * weight(23) = 20×23 -> 2 DSP
        fixed_type_weight cw = fixed_type_weight(c * w);
        fixed_type_weight sw = fixed_type_weight(s * w);
        // A += trig(20) * weight(23): 20×23 -> 2 DSP
        A_o += fixed_type_accum_ABC(c * cw);
        B_o += fixed_type_accum_ABC(s * sw);
        C_o += fixed_type_accum_ABC(c * sw);
        // D += radius(20) * weight(23): 20×23 -> 2 DSP
        D_o += fixed_type_accum_DE(r[i] * cw);
        E_o += fixed_type_accum_DE(r[i] * sw);
    }
}

// ── Stage 2b: A..E → hcx, hcy, rho ──────────────────────────────────────────
// Solver muls now use precisely-sized accum types -> 4 DSP each instead of 12
static void stage2b_solve(
    fixed_type_accum_ABC A_i, fixed_type_accum_ABC B_i, fixed_type_accum_ABC C_i,
    fixed_type_accum_DE  D_i, fixed_type_accum_DE  E_i,
    ap_uint<3>  nhits_i,
    fixed_type_hcxy &hcx_o, fixed_type_hcxy &hcy_o,
    fixed_type_hcxy &rho_o, ap_uint<3>  &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1

    nhits_o = nhits_i;
    if (nhits_i < 2) { hcx_o = hcy_o = rho_o = fixed_type_hcxy(0); return; }

    // A*B: accum_ABC(30) * accum_ABC(30) = 30×30 = 4 DSP (was 64×64 = 12)
    fixed_type_den den = fixed_type_den(FX_TWO) * fixed_type_den(
        fixed_type_den(A_i) * fixed_type_den(B_i) -
        fixed_type_den(C_i) * fixed_type_den(C_i));
    if (den == fixed_type_den(0)) { hcx_o = hcy_o = rho_o = fixed_type_hcxy(0); return; }

    // recip_den stored as ap_fixed<32,2>: keeps full 30 frac bits, range [-2,2] sufficient
    // Nx×recip: 64×32 = 6 DSP (was 64×64 = 12 DSP) -> saves 12 DSP
    fixed_type_recip_den recip_den = fixed_type_recip_den(
                                     call_den_recip_lut(fixed_type_general(den)));

    // Nx,Ny stored as ap_fixed<48,37>: covers 4.03e10, 48×32=4 DSP (vs 64×32=6 DSP)
    fixed_type_Nxy Nx = fixed_type_Nxy(D_i) * fixed_type_Nxy(B_i)
                      - fixed_type_Nxy(E_i) * fixed_type_Nxy(C_i);
    fixed_type_Nxy Ny = fixed_type_Nxy(E_i) * fixed_type_Nxy(A_i)
                      - fixed_type_Nxy(D_i) * fixed_type_Nxy(C_i);

    // hcx,hcy: Nx(48) × recip(32) = 4 DSP (was 6 DSP)
    fixed_type_hcxy hcx = fixed_type_hcxy(fixed_type_general(Nx) * fixed_type_general(recip_den));
    fixed_type_hcxy hcy = fixed_type_hcxy(fixed_type_general(Ny) * fixed_type_general(recip_den));

    // rho: hcxy(28) * hcxy(28) = 28×28 = 4 DSP (was 64×64 = 12)
    fixed_type_hcxy rho = fixed_type_hcxy(call_sqrt_lut(
        fixed_type_general(hcx * hcx) + fixed_type_general(hcy * hcy)));

    static const fixed_type_hcxy RHO_MAX = fixed_type_hcxy("3000.0");
    if (rho > RHO_MAX) { hcx_o = hcy_o = rho_o = fixed_type_hcxy(0); return; }

    hcx_o = hcx; hcy_o = hcy; rho_o = rho;
}

// ── Stage 2c_recip: compute recip(|hcx|) and recip(|hcy|) in parallel ─────────
// Separated from atan2 so each stage has a shorter critical path -> II=1
static void stage2c_recip(
    fixed_type_hcxy hcx_i, fixed_type_hcxy hcy_i,
    fixed_type_hcxy rho_i, ap_uint<3> nhits_i,
    fixed_type_general &recip_hcx_o, fixed_type_general &recip_hcy_o,
    fixed_type_hcxy &hcx_o, fixed_type_hcxy &hcy_o,
    fixed_type_hcxy &rho_o,  ap_uint<3> &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS ALLOCATION instances=call_small_recip_lut limit=2 function

    nhits_o = nhits_i; rho_o = rho_i;
    hcx_o = hcx_i; hcy_o = hcy_i;

    if (nhits_i < 2 || rho_i == fixed_type_hcxy(0)) {
        recip_hcx_o = recip_hcy_o = FX_ZERO; return;
    }

    fixed_type_general hcx    = fixed_type_general(hcx_i);
    fixed_type_general hcy    = fixed_type_general(hcy_i);
    fixed_type_general abs_hcx = (hcx < FX_ZERO) ? fixed_type_general(-hcx) : hcx;
    fixed_type_general abs_hcy = (hcy < FX_ZERO) ? fixed_type_general(-hcy) : hcy;

    recip_hcx_o = (abs_hcx != FX_ZERO) ? call_small_recip_lut(abs_hcx) : FX_ZERO;
    recip_hcy_o = (abs_hcy != FX_ZERO) ? call_small_recip_lut(abs_hcy) : FX_ZERO;
}

// ── Stage 2c_atan: atan2 lookup with pre-computed reciprocals ─────────────────
// Short critical path (atan_unit_lut only) -> II=1
static void stage2c_atan(
    fixed_type_general recip_hcx_i, fixed_type_general recip_hcy_i,
    fixed_type_hcxy hcx_i, fixed_type_hcxy hcy_i,
    fixed_type_hcxy rho_i, ap_uint<3> nhits_i,
    fixed_type_general &phi0_int_o,
    fixed_type_hcxy &rho_o, ap_uint<3> &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1

    nhits_o = nhits_i; rho_o = rho_i;

    if (nhits_i < 2 || rho_i == fixed_type_hcxy(0)) { phi0_int_o = FX_ZERO; return; }

    fixed_type_general hcx = fixed_type_general(hcx_i);
    fixed_type_general hcy = fixed_type_general(hcy_i);
    fixed_type_general phi0 = my_atan2_with_recip(hcy, hcx, recip_hcx_i, recip_hcy_i);
    if (phi0 < FX_ZERO) phi0 += FX_TWO_PI;
    phi0_int_o = phi0;
}

// ── Stage 3 ───────────────────────────────────────────────────────────────────
static void stage3_out(
    fixed_type_hcxy    rho_i, fixed_type_general phi0_int_i,
    fixed_type_general tsId[5], fixed_type_general useSL[5],
    ap_uint<3> nhits_i, ap_int<4> finderCharge,
    fixed_type_general &rho_o,  fixed_type_general &phi0_int_o,
    fixed_type_general &phi0_o, fixed_type_general &omega_o,
    fixed_type_general &chi2_o, ap_int<3> &chargeFit_o, ap_uint<3> &nhits_o
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=tsId  complete
#pragma HLS ARRAY_PARTITION variable=useSL complete

    nhits_o = nhits_i; chi2_o = FX_ZERO;

    if (nhits_i < 2 || rho_i == fixed_type_hcxy(0)) {
        rho_o = phi0_int_o = phi0_o = omega_o = FX_ZERO;
        chargeFit_o = 0; return;
    }

    rho_o = fixed_type_general(rho_i); phi0_int_o = phi0_int_i;

    ap_int<2> vote[5];
#pragma HLS ARRAY_PARTITION variable=vote complete
    for (int i = 0; i < 5; ++i) {
#pragma HLS UNROLL
        if (useSL[i] == FX_ZERO) { vote[i] = 0; continue; }
        // tsId is wire index ap_uint<10>, cast reduces mul width: 10×20 = 1 DSP
        fixed_type_step ts_n = (fixed_type_step)TS_PHI_STEP[i];
        ap_uint<10> tsId_int = (ap_uint<10>)tsId[i];
        fixed_type_general dPhi = fixed_type_general(tsId_int * ts_n) - phi0_int_i;
        if      (dPhi >  FX_PI)  dPhi -= FX_TWO_PI;
        else if (dPhi <= -FX_PI) dPhi += FX_TWO_PI;
        vote[i] = (dPhi > FX_ZERO) ? ap_int<2>(1)
                : (dPhi < FX_ZERO) ? ap_int<2>(-1) : ap_int<2>(0);
    }
    ap_int<4> sumCharge = vote[0]+vote[1]+vote[2]+vote[3]+vote[4];
    ap_int<3> chargeFit;
    if      (sumCharge > 0) chargeFit =  1;
    else if (sumCharge < 0) chargeFit = -1;
    else chargeFit = (finderCharge > 0) ? ap_int<3>(1)
                   : (finderCharge < 0) ? ap_int<3>(-1) : ap_int<3>(0);
    chargeFit_o = chargeFit;

    omega_o = (chargeFit == 0)
              ? FX_ZERO
              : fixed_type_general(fixed_type_general(chargeFit) * call_small_recip_lut(rho_o));
    phi0_o = wrap_to_npi_pi(phi0_int_i);
    (void)stored_track_phi0(phi0_int_i, chargeFit);
}

// ── Top-level kernel ──────────────────────────────────────────────────────────
void cdc2dfit_kernel(
    ap_uint<6> layer0, ap_uint<4> sl0_raw, ap_uint<10> wire0, ap_uint<2> lr0, ap_uint<10> priorityTime0, ap_uint<14> driftTick0,
    ap_uint<6> layer1, ap_uint<4> sl1_raw, ap_uint<10> wire1, ap_uint<2> lr1, ap_uint<10> priorityTime1, ap_uint<14> driftTick1,
    ap_uint<6> layer2, ap_uint<4> sl2_raw, ap_uint<10> wire2, ap_uint<2> lr2, ap_uint<10> priorityTime2, ap_uint<14> driftTick2,
    ap_uint<6> layer3, ap_uint<4> sl3_raw, ap_uint<10> wire3, ap_uint<2> lr3, ap_uint<10> priorityTime3, ap_uint<14> driftTick3,
    ap_uint<6> layer4, ap_uint<4> sl4_raw, ap_uint<10> wire4, ap_uint<2> lr4, ap_uint<10> priorityTime4, ap_uint<14> driftTick4,
    ap_int<16> T0, ap_int<4> finderCharge,
    fixed_type_general &rho,   fixed_type_general &phi0_internal,
    fixed_type_general &phi0,  fixed_type_general &omega,
    fixed_type_general &chi2,  ap_int<3> &chargeFit, ap_uint<3> &nhits
) {
#pragma HLS INTERFACE ap_none port=layer0
#pragma HLS INTERFACE ap_none port=sl0_raw
#pragma HLS INTERFACE ap_none port=wire0
#pragma HLS INTERFACE ap_none port=lr0
#pragma HLS INTERFACE ap_none port=priorityTime0
#pragma HLS INTERFACE ap_none port=driftTick0
#pragma HLS INTERFACE ap_none port=layer1
#pragma HLS INTERFACE ap_none port=sl1_raw
#pragma HLS INTERFACE ap_none port=wire1
#pragma HLS INTERFACE ap_none port=lr1
#pragma HLS INTERFACE ap_none port=priorityTime1
#pragma HLS INTERFACE ap_none port=driftTick1
#pragma HLS INTERFACE ap_none port=layer2
#pragma HLS INTERFACE ap_none port=sl2_raw
#pragma HLS INTERFACE ap_none port=wire2
#pragma HLS INTERFACE ap_none port=lr2
#pragma HLS INTERFACE ap_none port=priorityTime2
#pragma HLS INTERFACE ap_none port=driftTick2
#pragma HLS INTERFACE ap_none port=layer3
#pragma HLS INTERFACE ap_none port=sl3_raw
#pragma HLS INTERFACE ap_none port=wire3
#pragma HLS INTERFACE ap_none port=lr3
#pragma HLS INTERFACE ap_none port=priorityTime3
#pragma HLS INTERFACE ap_none port=driftTick3
#pragma HLS INTERFACE ap_none port=layer4
#pragma HLS INTERFACE ap_none port=sl4_raw
#pragma HLS INTERFACE ap_none port=wire4
#pragma HLS INTERFACE ap_none port=lr4
#pragma HLS INTERFACE ap_none port=priorityTime4
#pragma HLS INTERFACE ap_none port=driftTick4
#pragma HLS INTERFACE ap_none port=T0
#pragma HLS INTERFACE ap_none port=finderCharge
#pragma HLS INTERFACE ap_none port=rho
#pragma HLS INTERFACE ap_none port=phi0_internal
#pragma HLS INTERFACE ap_none port=phi0
#pragma HLS INTERFACE ap_none port=omega
#pragma HLS INTERFACE ap_none port=chi2
#pragma HLS INTERFACE ap_none port=chargeFit
#pragma HLS INTERFACE ap_none port=nhits
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS DATAFLOW

    ap_uint<6>  layer_a[5]  = { layer0,layer1,layer2,layer3,layer4 };
    ap_uint<4>  sl_raw_a[5] = { sl0_raw,sl1_raw,sl2_raw,sl3_raw,sl4_raw };
    ap_uint<10> wire_a[5]   = { wire0,wire1,wire2,wire3,wire4 };
    ap_uint<2>  lr_a[5]     = { lr0,lr1,lr2,lr3,lr4 };
    ap_uint<10> pTime_a[5]  = { priorityTime0,priorityTime1,priorityTime2,priorityTime3,priorityTime4 };
    ap_uint<14> dTick_a[5]  = { driftTick0,driftTick1,driftTick2,driftTick3,driftTick4 };
#pragma HLS ARRAY_PARTITION variable=layer_a  complete
#pragma HLS ARRAY_PARTITION variable=sl_raw_a complete
#pragma HLS ARRAY_PARTITION variable=wire_a   complete
#pragma HLS ARRAY_PARTITION variable=lr_a     complete
#pragma HLS ARRAY_PARTITION variable=pTime_a  complete
#pragma HLS ARRAY_PARTITION variable=dTick_a  complete

    fixed_type_general      phi_s[5], invErr_s[5], tsId_s[5], useSL_s[5];
    fixed_type_accum_ABC    A_s, B_s, C_s;
    fixed_type_accum_DE     D_s, E_s;
    fixed_type_hcxy         hcx_s, hcy_s, rho_s2b;           // stage2b output
    fixed_type_hcxy         hcx_s2cr, hcy_s2cr, rho_s2cr;   // stage2c_recip pass-through
    fixed_type_hcxy         rho_s2c;
    fixed_type_general      recip_hcx_s, recip_hcy_s;
    fixed_type_general      phi0_int_s;
    ap_uint<3>              nhits_s1, nhits_s2a, nhits_s2b, nhits_s2cr, nhits_s2c;
#pragma HLS ARRAY_PARTITION variable=phi_s    complete
#pragma HLS ARRAY_PARTITION variable=invErr_s complete
#pragma HLS ARRAY_PARTITION variable=tsId_s   complete
#pragma HLS ARRAY_PARTITION variable=useSL_s  complete

    stage1_hit   (layer_a, sl_raw_a, wire_a, lr_a, pTime_a, dTick_a, T0,
                  phi_s, invErr_s, tsId_s, useSL_s, nhits_s1);

    stage2a_accum(phi_s, invErr_s, nhits_s1,
                  A_s, B_s, C_s, D_s, E_s, nhits_s2a);

    stage2b_solve(A_s, B_s, C_s, D_s, E_s, nhits_s2a,
                  hcx_s, hcy_s, rho_s2b, nhits_s2b);

    stage2c_recip(hcx_s, hcy_s, rho_s2b, nhits_s2b,
                  recip_hcx_s, recip_hcy_s,
                  hcx_s2cr, hcy_s2cr, rho_s2cr, nhits_s2cr);

    stage2c_atan (recip_hcx_s, recip_hcy_s, hcx_s2cr, hcy_s2cr, rho_s2cr, nhits_s2cr,
                  phi0_int_s, rho_s2c, nhits_s2c);

    stage3_out   (rho_s2c, phi0_int_s, tsId_s, useSL_s, nhits_s2c, finderCharge,
                  rho, phi0_internal, phi0, omega, chi2, chargeFit, nhits);
}