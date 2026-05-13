#pragma once
// This header assumes fixed_type_general is already defined in cdc2dfit_kernel.h


static const fixed_type_general INV_RADIUS[5] = {
    fixed_type_general("0.050505050505"),  // 1/19.8
    fixed_type_general("0.024900398406"),  // 1/40.16
    fixed_type_general("0.016129032258"),  // 1/62.0
    fixed_type_general("0.011927480916"),  // 1/83.84
    fixed_type_general("0.009462528390")   // 1/105.68
};

static const fixed_type_general INV_WIRE_PHI_ERR[5] = {
    fixed_type_general("1174.998824054"),
    fixed_type_general("2509.976155217"),
    fixed_type_general("3875.842827249"),
    fixed_type_general("5239.991615657"),
    fixed_type_general("6605.019815060")
};

static const fixed_type_general INV_DRIFT_PHI_ERR[5] = {
    fixed_type_general("1174.998824054"),
    fixed_type_general("2509.976155217"),
    fixed_type_general("3875.842827249"),
    fixed_type_general("5239.991615657"),
    fixed_type_general("6605.019815060")
};


static const int DEN_MANT_LUT_SIZE = 256;

static const int DEN_EXP_MIN = -30;
static const int DEN_EXP_MAX =30;

// Mantissa LUT domain: m in [1, 2)
static const fixed_type_general DEN_MANT_LUT_X_MIN = fixed_type_general("1.000000000000000000e+00");
static const fixed_type_general DEN_MANT_LUT_X_MAX = fixed_type_general("2.000000000000000000e+00");
static const fixed_type_general DEN_MANT_LUT_STEP = fixed_type_general("3.906250000000000000e-03");
static const fixed_type_general DEN_MANT_LUT_STEP_INV = fixed_type_general("256.0");

static const fixed_type_general den_mant_recip_lut[DEN_MANT_LUT_SIZE] = {
    fixed_type_general("9.980506822612085216e-01"),
    fixed_type_general("9.941747572815533562e-01"),
    fixed_type_general("9.903288201160541648e-01"),
    fixed_type_general("9.865125240847784083e-01"),
    fixed_type_general("9.827255278310940145e-01"),
    fixed_type_general("9.789674952198852420e-01"),
    fixed_type_general("9.752380952380952372e-01"),
    fixed_type_general("9.715370018975332256e-01"),
    fixed_type_general("9.678638941398866047e-01"),
    fixed_type_general("9.642184557438794323e-01"),
    fixed_type_general("9.606003752345215752e-01"),
    fixed_type_general("9.570093457943925630e-01"),
    fixed_type_general("9.534450651769087459e-01"),
    fixed_type_general("9.499072356215213508e-01"),
    fixed_type_general("9.463955637707948121e-01"),
    fixed_type_general("9.429097605893186351e-01"),
    fixed_type_general("9.394495412844037219e-01"),
    fixed_type_general("9.360146252285191926e-01"),
    fixed_type_general("9.326047358834244072e-01"),
    fixed_type_general("9.292196007259527768e-01"),
    fixed_type_general("9.258589511754068413e-01"),
    fixed_type_general("9.225225225225225367e-01"),
    fixed_type_general("9.192100538599641268e-01"),
    fixed_type_general("9.159212880143112745e-01"),
    fixed_type_general("9.126559714795008382e-01"),
    fixed_type_general("9.094138543516874229e-01"),
    fixed_type_general("9.061946902654867131e-01"),
    fixed_type_general("9.029982363315696148e-01"),
    fixed_type_general("8.998242530755711233e-01"),
    fixed_type_general("8.966725043782837190e-01"),
    fixed_type_general("8.935427574171029841e-01"),
    fixed_type_general("8.904347826086956852e-01"),
    fixed_type_general("8.873483535528595700e-01"),
    fixed_type_general("8.842832469775474546e-01"),
    fixed_type_general("8.812392426850258476e-01"),
    fixed_type_general("8.782161234991423537e-01"),
    fixed_type_general("8.752136752136752129e-01"),
    fixed_type_general("8.722316865417376608e-01"),
    fixed_type_general("8.692699490662139095e-01"),
    fixed_type_general("8.663282571912013230e-01"),
    fixed_type_general("8.634064080944350295e-01"),
    fixed_type_general("8.605042016806723204e-01"),
    fixed_type_general("8.576214405360134130e-01"),
    fixed_type_general("8.547579298831385897e-01"),
    fixed_type_general("8.519134775374376245e-01"),
    fixed_type_general("8.490878938640132878e-01"),
    fixed_type_general("8.462809917355371692e-01"),
    fixed_type_general("8.434925864909390558e-01"),
    fixed_type_general("8.407224958949096605e-01"),
    fixed_type_general("8.379705400981997121e-01"),
    fixed_type_general("8.352365415986949815e-01"),
    fixed_type_general("8.325203252032520318e-01"),
    fixed_type_general("8.298217179902754870e-01"),
    fixed_type_general("8.271405492730210529e-01"),
    fixed_type_general("8.244766505636070830e-01"),
    fixed_type_general("8.218298555377206993e-01"),
    fixed_type_general("8.192000000000000393e-01"),
    fixed_type_general("8.165869218500797500e-01"),
    fixed_type_general("8.139904610492846304e-01"),
    fixed_type_general("8.114104595879556570e-01"),
    fixed_type_general("8.088467614533965122e-01"),
    fixed_type_general("8.062992125984251857e-01"),
    fixed_type_general("8.037676609105181003e-01"),
    fixed_type_general("8.012519561815336644e-01"),
    fixed_type_general("7.987519500780031478e-01"),
    fixed_type_general("7.962674961119751149e-01"),
    fixed_type_general("7.937984496124030898e-01"),
    fixed_type_general("7.913446676970633531e-01"),
    fixed_type_general("7.889060092449923234e-01"),
    fixed_type_general("7.864823348694316429e-01"),
    fixed_type_general("7.840735068912710881e-01"),
    fixed_type_general("7.816793893129770909e-01"),
    fixed_type_general("7.792998477929984347e-01"),
    fixed_type_general("7.769347496206373549e-01"),
    fixed_type_general("7.745839636913767201e-01"),
    fixed_type_general("7.722473604826546323e-01"),
    fixed_type_general("7.699248120300752340e-01"),
    fixed_type_general("7.676161919040479509e-01"),
    fixed_type_general("7.653213751868460646e-01"),
    fixed_type_general("7.630402384500745150e-01"),
    fixed_type_general("7.607726597325408235e-01"),
    fixed_type_general("7.585185185185184809e-01"),
    fixed_type_general("7.562776957163959146e-01"),
    fixed_type_general("7.540500736377024893e-01"),
    fixed_type_general("7.518355359765050983e-01"),
    fixed_type_general("7.496339677891654674e-01"),
    fixed_type_general("7.474452554744525079e-01"),
    fixed_type_general("7.452692867540029464e-01"),
    fixed_type_general("7.431059506531204617e-01"),
    fixed_type_general("7.409551374819102199e-01"),
    fixed_type_general("7.388167388167388161e-01"),
    fixed_type_general("7.366906474820144046e-01"),
    fixed_type_general("7.345767575322812437e-01"),
    fixed_type_general("7.324749642346208844e-01"),
    fixed_type_general("7.303851640513552290e-01"),
    fixed_type_general("7.283072546230441313e-01"),
    fixed_type_general("7.262411347517730986e-01"),
    fixed_type_general("7.241867043847242114e-01"),
    fixed_type_general("7.221438645980253757e-01"),
    fixed_type_general("7.201125175808720247e-01"),
    fixed_type_general("7.180925666199158286e-01"),
    fixed_type_general("7.160839160839160833e-01"),
    fixed_type_general("7.140864714086471166e-01"),
    fixed_type_general("7.121001390820583810e-01"),
    fixed_type_general("7.101248266296810163e-01"),
    fixed_type_general("7.081604426002766628e-01"),
    fixed_type_general("7.062068965517240837e-01"),
    fixed_type_general("7.042640990371389353e-01"),
    fixed_type_general("7.023319615912207992e-01"),
    fixed_type_general("7.004103967168262557e-01"),
    fixed_type_general("6.984993178717598949e-01"),
    fixed_type_general("6.965986394557822647e-01"),
    fixed_type_general("6.947082767978289830e-01"),
    fixed_type_general("6.928281461434371291e-01"),
    fixed_type_general("6.909581646423751389e-01"),
    fixed_type_general("6.890982503364737610e-01"),
    fixed_type_general("6.872483221476509696e-01"),
    fixed_type_general("6.854082998661311654e-01"),
    fixed_type_general("6.835781041388517831e-01"),
    fixed_type_general("6.817576564580559717e-01"),
    fixed_type_general("6.799468791500663523e-01"),
    fixed_type_general("6.781456953642384100e-01"),
    fixed_type_general("6.763540290620871920e-01"),
    fixed_type_general("6.745718050065876437e-01"),
    fixed_type_general("6.727989487516425893e-01"),
    fixed_type_general("6.710353866317169125e-01"),
    fixed_type_general("6.692810457516339406e-01"),
    fixed_type_general("6.675358539765319232e-01"),
    fixed_type_general("6.657997399219766077e-01"),
    fixed_type_general("6.640726329442282472e-01"),
    fixed_type_general("6.623544631306598207e-01"),
    fixed_type_general("6.606451612903225801e-01"),
    fixed_type_general("6.589446589446589231e-01"),
    fixed_type_general("6.572528883183568205e-01"),
    fixed_type_general("6.555697823303456850e-01"),
    fixed_type_general("6.538952745849297976e-01"),
    fixed_type_general("6.522292993630572910e-01"),
    fixed_type_general("6.505717916137230272e-01"),
    fixed_type_general("6.489226869455005930e-01"),
    fixed_type_general("6.472819216182048585e-01"),
    fixed_type_general("6.456494325346784358e-01"),
    fixed_type_general("6.440251572327043705e-01"),
    fixed_type_general("6.424090338770388486e-01"),
    fixed_type_general("6.408010012515644727e-01"),
    fixed_type_general("6.392009987515605562e-01"),
    fixed_type_general("6.376089663760896586e-01"),
    fixed_type_general("6.360248447204969180e-01"),
    fixed_type_general("6.344485749690210730e-01"),
    fixed_type_general("6.328800988875153966e-01"),
    fixed_type_general("6.313193588162762104e-01"),
    fixed_type_general("6.297662976629766485e-01"),
    fixed_type_general("6.282208588957055584e-01"),
    fixed_type_general("6.266829865361076557e-01"),
    fixed_type_general("6.251526251526251521e-01"),
    fixed_type_general("6.236297198538367503e-01"),
    fixed_type_general("6.221142162818954491e-01"),
    fixed_type_general("6.206060606060606055e-01"),
    fixed_type_general("6.191051995163240340e-01"),
    fixed_type_general("6.176115802171290303e-01"),
    fixed_type_general("6.161251504211793240e-01"),
    fixed_type_general("6.146458583433372924e-01"),
    fixed_type_general("6.131736526946107713e-01"),
    fixed_type_general("6.117084826762245742e-01"),
    fixed_type_general("6.102502979737782773e-01"),
    fixed_type_general("6.087990487514862714e-01"),
    fixed_type_general("6.073546856465006361e-01"),
    fixed_type_general("6.059171597633136175e-01"),
    fixed_type_general("6.044864226682408193e-01"),
    fixed_type_general("6.030624263839811094e-01"),
    fixed_type_general("6.016451233842537993e-01"),
    fixed_type_general("6.002344665885110953e-01"),
    fixed_type_general("5.988304093567251574e-01"),
    fixed_type_general("5.974329054842473230e-01"),
    fixed_type_general("5.960419091967403826e-01"),
    fixed_type_general("5.946573751451800227e-01"),
    fixed_type_general("5.932792584009269898e-01"),
    fixed_type_general("5.919075144508670894e-01"),
    fixed_type_general("5.905420991926182417e-01"),
    fixed_type_general("5.891829689298043737e-01"),
    fixed_type_general("5.878300803673938146e-01"),
    fixed_type_general("5.864833906071019731e-01"),
    fixed_type_general("5.851428571428571868e-01"),
    fixed_type_general("5.838084378563284105e-01"),
    fixed_type_general("5.824800910125141895e-01"),
    fixed_type_general("5.811577752553915843e-01"),
    fixed_type_general("5.798414496036240484e-01"),
    fixed_type_general("5.785310734463277038e-01"),
    fixed_type_general("5.772266065388951262e-01"),
    fixed_type_general("5.759280089988751961e-01"),
    fixed_type_general("5.746352413019080174e-01"),
    fixed_type_general("5.733482642777155691e-01"),
    fixed_type_general("5.720670391061452031e-01"),
    fixed_type_general("5.707915273132664336e-01"),
    fixed_type_general("5.695216907675194618e-01"),
    fixed_type_general("5.682574916759156602e-01"),
    fixed_type_general("5.669988925802879054e-01"),
    fixed_type_general("5.657458563535912033e-01"),
    fixed_type_general("5.644983461962513882e-01"),
    fixed_type_general("5.632563256325632262e-01"),
    fixed_type_general("5.620197585071350366e-01"),
    fixed_type_general("5.607886089813800545e-01"),
    fixed_type_general("5.595628415300546443e-01"),
    fixed_type_general("5.583424209378408110e-01"),
    fixed_type_general("5.571273122959738977e-01"),
    fixed_type_general("5.559174809989142485e-01"),
    fixed_type_general("5.547128927410617250e-01"),
    fixed_type_general("5.535135135135135220e-01"),
    fixed_type_general("5.523193096008629510e-01"),
    fixed_type_general("5.511302475780408550e-01"),
    fixed_type_general("5.499462943071965482e-01"),
    fixed_type_general("5.487674169346195008e-01"),
    fixed_type_general("5.475935828877005473e-01"),
    fixed_type_general("5.464247598719317311e-01"),
    fixed_type_general("5.452609158679446733e-01"),
    fixed_type_general("5.441020191285865781e-01"),
    fixed_type_general("5.429480381760339869e-01"),
    fixed_type_general("5.417989417989418355e-01"),
    fixed_type_general("5.406546990496303717e-01"),
    fixed_type_general("5.395152792413066001e-01"),
    fixed_type_general("5.383806519453206985e-01"),
    fixed_type_general("5.372507869884575182e-01"),
    fixed_type_general("5.361256544502618349e-01"),
    fixed_type_general("5.350052246603971273e-01"),
    fixed_type_general("5.338894681960375532e-01"),
    fixed_type_general("5.327783558792924534e-01"),
    fixed_type_general("5.316718587746624980e-01"),
    fixed_type_general("5.305699481865284728e-01"),
    fixed_type_general("5.294725956566701530e-01"),
    fixed_type_general("5.283797729618162631e-01"),
    fixed_type_general("5.272914521112255226e-01"),
    fixed_type_general("5.262076053442960033e-01"),
    fixed_type_general("5.251282051282051277e-01"),
    fixed_type_general("5.240532241555783122e-01"),
    fixed_type_general("5.229826353421859197e-01"),
    fixed_type_general("5.219164118246687467e-01"),
    fixed_type_general("5.208545269582909309e-01"),
    fixed_type_general("5.197969543147208382e-01"),
    fixed_type_general("5.187436676798379276e-01"),
    fixed_type_general("5.176946410515672614e-01"),
    fixed_type_general("5.166498486377396615e-01"),
    fixed_type_general("5.156092648539778445e-01"),
    fixed_type_general("5.145728643216080922e-01"),
    fixed_type_general("5.135406218655967914e-01"),
    fixed_type_general("5.125125125125125081e-01"),
    fixed_type_general("5.114885114885114881e-01"),
    fixed_type_general("5.104685942173479152e-01"),
    fixed_type_general("5.094527363184079283e-01"),
    fixed_type_general("5.084409136047666200e-01"),
    fixed_type_general("5.074331020812685722e-01"),
    fixed_type_general("5.064292779426310398e-01"),
    fixed_type_general("5.054294175715695614e-01"),
    fixed_type_general("5.044334975369457741e-01"),
    fixed_type_general("5.034414945919370998e-01"),
    fixed_type_general("5.024533856722276370e-01"),
    fixed_type_general("5.014691478942213676e-01"),
    fixed_type_general("5.004887585532746819e-01")
};

static const int DEN_EXP_LUT_SIZE = DEN_EXP_MAX - DEN_EXP_MIN + 1;

static const fixed_type_general den_exp_recip_lut[DEN_EXP_LUT_SIZE] = {
    fixed_type_general("1073741824.0"),   // 2^30
    fixed_type_general("536870912.0"),    // 2^29
    fixed_type_general("268435456.0"),    // 2^28
    fixed_type_general("134217728.0"),    // 2^27
    fixed_type_general("67108864.0"),     // 2^26
    fixed_type_general("33554432.0"),     // 2^25
    fixed_type_general("16777216.0"),     // 2^24
    fixed_type_general("8388608.0"),      // 2^23
    fixed_type_general("4194304.0"),      // 2^22
    fixed_type_general("2097152.0"),      // 2^21
    fixed_type_general("1048576.0"),      // 2^20
    fixed_type_general("524288.0"),       // 2^19
    fixed_type_general("262144.0"),       // 2^18
    fixed_type_general("131072.0"),       // 2^17
    fixed_type_general("65536.0"),        // 2^16
    fixed_type_general("32768.0"),        // 2^15
    fixed_type_general("16384.0"),        // 2^14
    fixed_type_general("8192.0"),         // 2^13
    fixed_type_general("4096.0"),         // 2^12
    fixed_type_general("2048.0"),         // 2^11
    fixed_type_general("1024.0"),         // 2^10
    fixed_type_general("512.0"),          // 2^9
    fixed_type_general("256.0"),          // 2^8
    fixed_type_general("128.0"),          // 2^7
    fixed_type_general("64.0"),           // 2^6
    fixed_type_general("32.0"),           // 2^5
    fixed_type_general("16.0"),           // 2^4
    fixed_type_general("8.0"),            // 2^3
    fixed_type_general("4.0"),            // 2^2
    fixed_type_general("2.0"),            // 2^1
    fixed_type_general("1.0"),            // 2^0
    fixed_type_general("0.5"),            // 2^-1
    fixed_type_general("0.25"),           // 2^-2
    fixed_type_general("0.125"),          // 2^-3
    fixed_type_general("0.0625"),         // 2^-4
    fixed_type_general("0.03125"),        // 2^-5
    fixed_type_general("0.015625"),       // 2^-6
    fixed_type_general("0.0078125"),      // 2^-7
    fixed_type_general("0.00390625"),     // 2^-8
    fixed_type_general("0.001953125"),    // 2^-9
    fixed_type_general("0.0009765625"),   // 2^-10
    fixed_type_general("0.00048828125"),  // 2^-11
    fixed_type_general("0.000244140625"), // 2^-12
    fixed_type_general("0.0001220703125"),// 2^-13
    fixed_type_general("0.00006103515625"), // 2^-14
    fixed_type_general("0.000030517578125"), // 2^-15
    fixed_type_general("0.0000152587890625"), // 2^-16
    fixed_type_general("7.62939453125e-6"),   // 2^-17
    fixed_type_general("3.814697265625e-6"),  // 2^-18
    fixed_type_general("1.9073486328125e-6"), // 2^-19
    fixed_type_general("9.5367431640625e-7"), // 2^-20
    fixed_type_general("4.76837158203125e-7"), // 2^-21
    fixed_type_general("2.384185791015625e-7"), // 2^-22
    fixed_type_general("1.1920928955078125e-7"), // 2^-23
    fixed_type_general("5.9604644775390625e-8"), // 2^-24
    fixed_type_general("2.98023223876953125e-8"), // 2^-25
    fixed_type_general("1.490116119384765625e-8"), // 2^-26
    fixed_type_general("7.450580596923828125e-9"), // 2^-27
    fixed_type_general("3.7252902984619140625e-9"), // 2^-28
    fixed_type_general("1.86264514923095703125e-9"), // 2^-29
    fixed_type_general("9.31322574615478515625e-10") // 2^-30
};

static inline fixed_type_general call_den_recip_lut(fixed_type_general den) {
#pragma HLS INLINE
#pragma HLS BIND_STORAGE variable=den_mant_recip_lut type=ROM_NP impl=LUTRAM
#pragma HLS BIND_STORAGE variable=den_exp_recip_lut  type=ROM_2P impl=LUTRAM

    fixed_type_general abs_den = den;
    if (abs_den < FX_ZERO) abs_den = -abs_den;
    if (abs_den == FX_ZERO) return FX_ZERO;

    fixed_type_general m = abs_den;
    int exp_shift = 0;

    for (int k = 0; k < 31; ++k) {
#pragma HLS UNROLL
        if ((m >= FX_TWO) && (exp_shift < DEN_EXP_MAX)) {
            m *= FX_HALF;
            exp_shift++;
        }
    }

    for (int k = 0; k < 31; ++k) {
#pragma HLS UNROLL
        if ((m < FX_ONE) && (exp_shift > DEN_EXP_MIN)) {
            m *= FX_TWO;
            exp_shift--;
        }
    }

    int midx = (int)((m - DEN_MANT_LUT_X_MIN) * DEN_MANT_LUT_STEP_INV + FX_HALF);
    if (midx < 0) midx = 0;
    if (midx >= DEN_MANT_LUT_SIZE) midx = DEN_MANT_LUT_SIZE - 1;

    int eidx = exp_shift - DEN_EXP_MIN;
    if (eidx < 0) eidx = 0;
    if (eidx >= DEN_EXP_LUT_SIZE) eidx = DEN_EXP_LUT_SIZE - 1;

    // 32×32 = 4 DSP (was 64×64 = 12 DSP); range [0.5,2]×[1/2^30,2^30] -> output small
    typedef ap_fixed<32,2> recip_dn_t;
    recip_dn_t recip_n = recip_dn_t(den_mant_recip_lut[midx]) * recip_dn_t(den_exp_recip_lut[eidx]);
    fixed_type_general recip = fixed_type_general(recip_n);
    if (den < FX_ZERO) recip = -recip;
    return recip;
}
// Reciprocal LUT for small inputs (|hcx|, |hcy|, rho): range [0, 3000]
// Uses 12 normalization iterations instead of 31, saving ~380 LUT/FF.
static const int SMALL_EXP_MIN = -12;
static const int SMALL_EXP_MAX =  12;
static const int SMALL_EXP_LUT_SIZE = SMALL_EXP_MAX - SMALL_EXP_MIN + 1;

static const fixed_type_general small_exp_recip_lut[SMALL_EXP_LUT_SIZE] = {
    fixed_type_general("4.096000000000000000e+03"),  // 2^(--12)
    fixed_type_general("2.048000000000000000e+03"),  // 2^(--11)
    fixed_type_general("1.024000000000000000e+03"),  // 2^(--10)
    fixed_type_general("5.120000000000000000e+02"),  // 2^(--9)
    fixed_type_general("2.560000000000000000e+02"),  // 2^(--8)
    fixed_type_general("1.280000000000000000e+02"),  // 2^(--7)
    fixed_type_general("6.400000000000000000e+01"),  // 2^(--6)
    fixed_type_general("3.200000000000000000e+01"),  // 2^(--5)
    fixed_type_general("1.600000000000000000e+01"),  // 2^(--4)
    fixed_type_general("8.000000000000000000e+00"),  // 2^(--3)
    fixed_type_general("4.000000000000000000e+00"),  // 2^(--2)
    fixed_type_general("2.000000000000000000e+00"),  // 2^(--1)
    fixed_type_general("1.000000000000000000e+00"),  // 2^(-0)
    fixed_type_general("5.000000000000000000e-01"),  // 2^(-1)
    fixed_type_general("2.500000000000000000e-01"),  // 2^(-2)
    fixed_type_general("1.250000000000000000e-01"),  // 2^(-3)
    fixed_type_general("6.250000000000000000e-02"),  // 2^(-4)
    fixed_type_general("3.125000000000000000e-02"),  // 2^(-5)
    fixed_type_general("1.562500000000000000e-02"),  // 2^(-6)
    fixed_type_general("7.812500000000000000e-03"),  // 2^(-7)
    fixed_type_general("3.906250000000000000e-03"),  // 2^(-8)
    fixed_type_general("1.953125000000000000e-03"),  // 2^(-9)
    fixed_type_general("9.765625000000000000e-04"),  // 2^(-10)
    fixed_type_general("4.882812500000000000e-04"),  // 2^(-11)
    fixed_type_general("2.441406250000000000e-04"),  // 2^(-12)
};

static inline fixed_type_general call_small_recip_lut(fixed_type_general x) {
#pragma HLS INLINE
#pragma HLS BIND_STORAGE variable=den_mant_recip_lut  type=ROM_NP impl=LUTRAM
#pragma HLS BIND_STORAGE variable=small_exp_recip_lut type=ROM_NP impl=LUTRAM

    if (x == FX_ZERO) return FX_ZERO;

    fixed_type_general m = (x < FX_ZERO) ? fixed_type_general(-x) : x;

    // Priority encoder: O(1) combinational mux tree instead of sequential loop.
    // Range [0, 3000] -> exp_shift in [-2, 11].
    int exp_shift;
    if      (m >= fixed_type_general("2048.0")) exp_shift = 11;
    else if (m >= fixed_type_general("1024.0")) exp_shift = 10;
    else if (m >= fixed_type_general( "512.0")) exp_shift =  9;
    else if (m >= fixed_type_general( "256.0")) exp_shift =  8;
    else if (m >= fixed_type_general( "128.0")) exp_shift =  7;
    else if (m >= fixed_type_general(  "64.0")) exp_shift =  6;
    else if (m >= fixed_type_general(  "32.0")) exp_shift =  5;
    else if (m >= fixed_type_general(  "16.0")) exp_shift =  4;
    else if (m >= fixed_type_general(   "8.0")) exp_shift =  3;
    else if (m >= fixed_type_general(   "4.0")) exp_shift =  2;
    else if (m >= fixed_type_general(   "2.0")) exp_shift =  1;
    else if (m >= fixed_type_general(   "1.0")) exp_shift =  0;
    else if (m >= fixed_type_general(   "0.5")) exp_shift = -1;
    else                                        exp_shift = -2;

    int eidx = exp_shift - SMALL_EXP_MIN;
    if (eidx < 0)                    eidx = 0;
    if (eidx >= SMALL_EXP_LUT_SIZE)  eidx = SMALL_EXP_LUT_SIZE - 1;

    // mant = m * 2^(-exp_shift) in [1, 2)
    fixed_type_general mant = m * small_exp_recip_lut[eidx];

    int midx = (int)((mant - DEN_MANT_LUT_X_MIN) * DEN_MANT_LUT_STEP_INV);
    if (midx < 0)                  midx = 0;
    if (midx >= DEN_MANT_LUT_SIZE) midx = DEN_MANT_LUT_SIZE - 1;

    // narrow 18x18 mul -> 1 DSP
    typedef ap_fixed<18,2> rn_t;
    rn_t recip_n = rn_t(den_mant_recip_lut[midx]) * rn_t(small_exp_recip_lut[eidx]);
    fixed_type_general recip = fixed_type_general(recip_n);

    if (x < FX_ZERO) recip = -recip;
    return recip;
}