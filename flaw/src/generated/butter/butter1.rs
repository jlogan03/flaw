//! Butterworth filter of order 1.
//! Region of validity: cutoff ratio from 1.00e-04 to 4.00e-01 .
//! This file is autogenerated.
#![allow(clippy::style)]

use crate::SisoIirFilter;

/// Minimum tabulated cutoff ratio
#[allow(dead_code)]
pub const MIN_CUTOFF_RATIO: f64 = 0.0001;

/// Maximum tabulated cutoff ratio
#[allow(dead_code)]
pub const MAX_CUTOFF_RATIO: f64 = 0.4;

/// Initialise a Butterworth filter of order 1 by interpolating the coefficients from stored tables.
/// Cutoff ratio is the dimensionless ratio of the cutoff frequency to the sampling frequency.
/// Region of validity: cutoff ratio from 1.00e-04 to 4.00e-01
pub fn butter1(cutoff_ratio: f64) -> Result<SisoIirFilter<1>, &'static str> {
    let avals = &[&AVALS[0][..]];
    let cvals = &[&CVALS[0][..]];
    SisoIirFilter::new_interpolated(cutoff_ratio, &LOG10_CUTOFF_RATIOS, avals, cvals, &DVALS)
}

/// [dimensionless] Log base-10 of cutoff ratios, to improve float precision during interpolation
#[rustfmt::skip]
static LOG10_CUTOFF_RATIOS: [f64; 100] = [-4.0, -3.963615555643152, -3.9272311112863036, -3.8908466669294555, -3.8544622225726077, -3.8180777782157596, -3.7816933338589114, -3.745308889502063, -3.708924445145215, -3.672540000788367, -3.636155556431519, -3.599771112074671, -3.5633866677178228, -3.5270022233609746, -3.4906177790041264, -3.454233334647278, -3.41784889029043, -3.3814644459335823, -3.345080001576734, -3.308695557219886, -3.2723111128630378, -3.2359266685061896, -3.199542224149342, -3.1631577797924937, -3.1267733354356455, -3.0903888910787973, -3.054004446721949, -3.017620002365101, -2.981235558008253, -2.9448511136514046, -2.9084666692945564, -2.8720822249377087, -2.8356977805808605, -2.7993133362240123, -2.7629288918671646, -2.7265444475103164, -2.6901600031534683, -2.65377555879662, -2.617391114439772, -2.5810066700829237, -2.5446222257260755, -2.5082377813692274, -2.471853337012379, -2.435468892655531, -2.3990844482986833, -2.362700003941835, -2.326315559584987, -2.289931115228139, -2.253546670871291, -2.217162226514443, -2.1807777821575947, -2.1443933378007465, -2.1080088934438983, -2.07162444908705, -2.035240004730202, -1.9988555603733542, -1.962471116016506, -1.9260866716596579, -1.8897022273028097, -1.8533177829459615, -1.8169333385891133, -1.7805488942322656, -1.7441644498754174, -1.7077800055185692, -1.671395561161721, -1.6350111168048729, -1.5986266724480247, -1.5622422280911765, -1.5258577837343288, -1.4894733393774806, -1.4530888950206324, -1.4167044506637843, -1.380320006306936, -1.343935561950088, -1.3075511175932402, -1.271166673236392, -1.2347822288795438, -1.1983977845226956, -1.1620133401658475, -1.1256288958089993, -1.089244451452151, -1.0528600070953034, -1.0164755627384552, -0.980091118381607, -0.9437066740247588, -0.9073222296679107, -0.8709377853110625, -0.8345533409542147, -0.7981688965973666, -0.7617844522405184, -0.7254000078836702, -0.689015563526822, -0.6526311191699739, -0.6162466748131261, -0.5798622304562779, -0.5434777860994298, -0.5070933417425816, -0.4707088973857334, -0.4343244530288852, -0.3979400086720376];

/// State-Space `D` 1x1 matrix
#[rustfmt::skip]
static DVALS: [f64; 100] = [0.0003140606106404424, 0.0003414962968963103, 0.0003713278181914204, 0.00040376422981415743, 0.0004390327991953734, 0.0004773805870978057, 0.0005190761651035754, 0.0005644114809694962, 0.0006137038843692006, 0.0006672983265622219, 0.0007255697486273807, 0.0007889256740754595, 0.0008578090229184327, 0.0009327011656236573, 0.0010141252368252938, 0.0011026497302055044, 0.0011988923975977627, 0.0013035244771064777, 0.0014172752768826736, 0.0015409371431453142, 0.001675370843091152, 0.0018215113954900622, 0.0019803743840129924, 0.0021530627906783525, 0.002340774389219106, 0.0025448097406523397, 0.002766580835855939, 0.00300762043249806, 0.003269592136192098, 0.003554301278222749, 0.0038637066445576114, 0.004199933113063084, 0.00456528525780901, 0.004962261980985006, 0.005393572234156909, 0.0058621518912382745, 0.006371181835491648, 0.006924107321934546, 0.007524658674504007, 0.008176873373999836, 0.008885119587914134, 0.009654121186461745, 0.010488984280112825, 0.011395225302316065, 0.012378800646471403, 0.013446137848110542, 0.014604168281186089, 0.015860361310847827, 0.017222759813577338, 0.018700016938550355, 0.020301433941125848, 0.022036998870017057, 0.02391742583370028, 0.02595419450885875, 0.028159589484299793, 0.030546738958338787, 0.033129652227112276, 0.03592325531728484, 0.03894342403157888, 0.04220701359290823, 0.04573188399730397, 0.04953692012344807, 0.053642045605420616, 0.05806822946524479, 0.06283748453539818, 0.06797285679374275, 0.07349840490237924, 0.07943916950901325, 0.08582113225917866, 0.09267116500826003, 0.10001697044544343, 0.1078870162829709, 0.11631046636273983, 0.12531711353219255, 0.1349373209919857, 0.14520198107640267, 0.15614250316325787, 0.1677908457118818, 0.1801796114157945, 0.1933422293017646, 0.20731325356110641, 0.2221288163426111, 0.23782728125089758, 0.2544501567744438, 0.2720433456989697, 0.29065882989672914, 0.3103569230893974, 0.3312092725856692, 0.3533028631221499, 0.37674538568335225, 0.4016725046652589, 0.42825782740709384, 0.45672682010890364, 0.48737664870816877, 0.5206051874182946, 0.556954690524724, 0.597179804000546, 0.6423577279319551, 0.6940750722389863, 0.7547627247472144];

/// State-Space `A` matrix, first row
#[rustfmt::skip]
static AVALS: [[f64; 100]; 1] = [[0.9993718787787192, 0.9993170074062072, 0.9992573443636171, 0.9991924715403716, 0.9991219344016092, 0.9990452388258045, 0.9989618476697929, 0.9988711770380609, 0.9987725922312616, 0.9986654033468756, 0.9985488605027453, 0.998422148651849, 0.998284381954163, 0.9981345976687526, 0.9979717495263493, 0.9977947005395891, 0.9976022152048044, 0.9973929510457872, 0.9971654494462346, 0.9969181257137093, 0.9966492583138177, 0.9963569772090198, 0.996039251231974, 0.9956938744186431, 0.9953184512215618, 0.9949103805186955, 0.9944668383282881, 0.9939847591350041, 0.9934608157276158, 0.9928913974435546, 0.9922725867108848, 0.9916001337738737, 0.990869429484382, 0.9900754760380301, 0.9892128555316863, 0.9882756962175234, 0.9872576363290169, 0.986151785356131, 0.984950682650992, 0.9836462532520002, 0.9822297608241716, 0.9806917576270764, 0.9790220314397744, 0.9772095493953679, 0.9752423987070571, 0.973107724303779, 0.9707916634376279, 0.9682792773783044, 0.9655544803728453, 0.9625999661228992, 0.9593971321177482, 0.9559260022599659, 0.9521651483325996, 0.9480916109822825, 0.9436808210314003, 0.9389065220833223, 0.9337406955457753, 0.9281534893654304, 0.9221131519368422, 0.9155859728141835, 0.908536232005392, 0.900926159753104, 0.8927159087891587, 0.8838635410695104, 0.8743250309292037, 0.8640542864125145, 0.8530031901952415, 0.8411216609819735, 0.8283577354816427, 0.8146576699834798, 0.7999660591091131, 0.7842259674340581, 0.7673790672745204, 0.7493657729356149, 0.7301253580160285, 0.7095960378471948, 0.6877149936734842, 0.6644183085762364, 0.6396407771684112, 0.6133155413964708, 0.5853734928777872, 0.5557423673147779, 0.5243454374982048, 0.49109968645111235, 0.4559133086020605, 0.41868234020654177, 0.3792861538212053, 0.3375814548286617, 0.2933942737557001, 0.24650922863329552, 0.19665499066948214, 0.14348434518581235, 0.08654635978219286, 0.025246702583662326, -0.04121037483658935, -0.11390938104944813, -0.19435960800109217, -0.28471545586391017, -0.38815014447797264, -0.5095254494944288]];

/// State-Space `C` vector
#[rustfmt::skip]
static CVALS: [[f64; 100]; 1] = [[0.0006279239531465732, 0.0006827593543510328, 0.0007423798676857152, 0.00080720240852176, 0.0008776800987932081, 0.0009543053897457358, 0.0010376134500767935, 0.001128185841299292, 0.0012266545038230217, 0.0013337060790111783, 0.001450086594334515, 0.0015766065407124882, 0.001714146373197265, 0.001863662468318603, 0.0020261935736586552, 0.0022028677875559643, 0.0023949101092334895, 0.002603650602088124, 0.0028305332153444216, 0.0030771253117323786, 0.003345127951258544, 0.003636386983452324, 0.003952905002624275, 0.004296854222595497, 0.004670590328955764, 0.005076667368072442, 0.005517853732669227, 0.00599714930366416, 0.006517803806910097, 0.007083336441292746, 0.007697556831044826, 0.008364587349817761, 0.009088886856647684, 0.009875275874034153, 0.010728963225423683, 0.011655574132884653, 0.0126611797550215, 0.013752328119457757, 0.014936076372672438, 0.016220024231650816, 0.017612348475645197, 0.01912183826115771, 0.020757930977768742, 0.02253074828525104, 0.024451131882052644, 0.026530678450160303, 0.028781773100001774, 0.03121762049987438, 0.03385227271596233, 0.03670065261009657, 0.0397785714421199, 0.04310273910163985, 0.046690765150379454, 0.05056114859251015, 0.05473325400875101, 0.059227271394699785, 0.06406415674084573, 0.069265550089388, 0.07485366751255104, 0.08085116319295257, 0.0872809575667222, 0.09416602733626239, 0.10152915309737318, 0.10939262038403294, 0.11777787014532357, 0.1267050950660802, 0.13619277875837033, 0.14625717571346303, 0.15691173103386247, 0.16816644036854372, 0.18002715213671744, 0.19249481600105772, 0.2055646835544436, 0.21922546917630423, 0.233458480790983, 0.24823673153578135, 0.2635240437383397, 0.27927415561434654, 0.29542983809169565, 0.31192202334077695, 0.3286689369180296, 0.3455752105856833, 0.3625309310874081, 0.37941054898380955, 0.3960715275197611, 0.41235254899958684, 0.4280710067597585, 0.44301938067788216, 0.4569599000636826, 0.46961660009930883, 0.4806634073223929, 0.4897061213432993, 0.4962548638042257, 0.499681302004326, 0.49915085250291386, 0.49351232645446574, 0.4811121713888309, 0.4594685545961029, 0.4246697326708645, 0.37019190815875014]];

#[cfg(feature = "std")]
#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    static CUTOFF_TEST_INPUT: [f32; 101] = [0.0, 0.5877852522924732, -0.9510565162951536, 0.9510565162951535, -0.5877852522924728, -4.898587196589413e-16, 0.5877852522924736, -0.9510565162951538, 0.9510565162951533, -0.5877852522924725, -9.797174393178826e-16, 0.587785252292474, -0.951056516295154, 0.9510565162951532, -0.587785252292472, -1.4695761589768238e-15, 0.5877852522924744, -0.9510565162951541, 0.951056516295153, -0.5877852522924716, -1.959434878635765e-15, 0.5877852522924748, -0.9510565162951542, 0.9510565162951529, -0.5877852522924712, -2.4492935982947065e-15, 0.5877852522924751, -0.9510565162951544, 0.9510565162951528, -0.5877852522924708, -2.9391523179536475e-15, 0.5877852522924756, -0.9510565162951545, 0.9510565162951525, -0.5877852522924705, -3.429011037612589e-15, 0.587785252292476, -0.9510565162951546, 0.9510565162951524, -0.58778525229247, -3.91886975727153e-15, 0.5877852522924764, -0.9510565162951549, 0.9510565162951523, -0.5877852522924697, -4.408728476930472e-15, 0.5877852522924768, -0.951056516295155, 0.9510565162951521, -0.5877852522924693, -4.898587196589413e-15, 0.5877852522924887, -0.9510565162951552, 0.9510565162951563, -0.5877852522924688, -1.9599300631450357e-14, 0.5877852522924776, -0.9510565162951509, 0.9510565162951519, -0.5877852522924569, -5.878304635907295e-15, 0.5877852522924665, -0.9510565162951554, 0.9510565162951473, -0.587785252292468, 7.842691359635767e-15, 0.5877852522924784, -0.95105651629516, 0.9510565162951515, -0.5877852522924791, -6.858022075225178e-15, 0.5877852522924902, -0.9510565162951558, 0.9510565162951558, -0.5877852522924673, -2.1558735510086122e-14, 0.5877852522924791, -0.9510565162951515, 0.9510565162951512, -0.5877852522924554, -7.83773951454306e-15, 0.587785252292468, -0.9510565162951561, 0.9510565162951466, -0.5877852522924665, 5.883256481000002e-15, 0.5877852522924799, -0.9510565162951606, 0.9510565162951509, -0.5877852522924776, -8.817456953860943e-15, 0.5877852522924918, -0.9510565162951563, 0.9510565162951552, -0.5877852522924657, -2.3518170388721888e-14, 0.5877852522924807, -0.9510565162951521, 0.9510565162951506, -0.5877852522924538, -9.797174393178826e-15];
    static CUTOFF_TEST_OUTPUT: [f32; 101] = [0.0, 0.44363839858649595, -0.5002286634938007, 0.25487923461667605, 0.14431615245617008, -0.5171711510360325, 0.7071502617836815, -0.6344946640364608, 0.323291178894994, 0.1094585257969724, -0.49941030314419593, 0.6981006577781912, -0.6298836604878166, 0.3209417552392504, 0.11065561694121828, -0.5000202515475538, 0.6984114420125804, -0.6300420129645394, 0.3210224398561308, 0.11061450607553514, -0.4999993045152374, 0.6984007689665237, -0.6300365747759498, 0.3210196689606451, 0.1106159179173033, -0.5000000238845491, 0.6984011355034956, -0.6300367615358649, 0.3210197641195746, 0.11061586943140711, -0.49999999917975113, 0.6984011229157723, -0.6300367551220993, 0.32101976085159745, 0.11061587109652476, -0.500000000028171, 0.6984011233480638, -0.6300367553423626, 0.3210197609638271, 0.11061587103934117, -0.49999999999903455, 0.6984011233332178, -0.6300367553347981, 0.3210197609599726, 0.11061587104130527, -0.5000000000000354, 0.6984011233337277, -0.6300367553350577, 0.32101976096010454, 0.1106158710412381, -0.5000000000000013, 0.6984011233337191, -0.6300367553350443, 0.32101976096010076, 0.11061587104124365, -0.5000000000000149, 0.6984011233337064, -0.6300367553350431, 0.3210197609601001, 0.11061587104124965, -0.4999999999999985, 0.6984011233337, -0.6300367553350517, 0.32101976096009754, 0.11061587104123904, -0.49999999999999123, 0.6984011233337157, -0.6300367553350541, 0.32101976096009865, 0.11061587104123344, -0.5000000000000079, 0.698401123333722, -0.6300367553350451, 0.3210197609601003, 0.1106158710412446, -0.5000000000000157, 0.6984011233337064, -0.6300367553350426, 0.3210197609600987, 0.11061587104125092, -0.4999999999999995, 0.6984011233337002, -0.6300367553350512, 0.32101976096009627, 0.11061587104124032, -0.4999999999999922, 0.6984011233337158, -0.6300367553350534, 0.32101976096009727, 0.11061587104123471, -0.5000000000000088, 0.6984011233337222, -0.6300367553350446, 0.32101976096009915, 0.11061587104124593, -0.5000000000000165, 0.6984011233337066, -0.6300367553350418, 0.32101976096009766, 0.11061587104125231, -0.5000000000000006];
    const STEP_TEST_MIN_OUTPUT: f32 = 1.0000000000000966;
    const STEP_TEST_MAX_OUTPUT: f32 = 1.1249546329059645;

    #[test]
    fn test() {
        let order = 1;
        println!("order {order}");
        let mut filter = butter1(0.4).unwrap();
        let out = (0..CUTOFF_TEST_INPUT.len()).map(|i| {filter.update(CUTOFF_TEST_INPUT[i])}).collect::<Vec<f32>>();
        // Check overall match to reference output to catch phase error, etc
        (0..CUTOFF_TEST_INPUT.len()).for_each(|i| { let expected = CUTOFF_TEST_OUTPUT[i]; let rel_err = (out[i] - expected).abs() / expected.abs().max(1e-4); assert!(rel_err < 0.05); });
        // Check approximate attenuation at cutoff frequency; should be -3dB or 1/sqrt(2) magnitude
        let maxmag = out.iter().fold(0.0_f32, |a, b| a.abs().max(b.abs()));
        let attenuation_rel_err = (maxmag - 0.707).abs() / 0.707;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);
        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter1(MIN_CUTOFF_RATIO).unwrap();
        (0..99999).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        let step_min_rel_err = (step_min_final - STEP_TEST_MIN_OUTPUT).abs() / STEP_TEST_MIN_OUTPUT;
        println!("order {order} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);
        let mut filtermax = butter1(MAX_CUTOFF_RATIO).unwrap();
        (0..1).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - STEP_TEST_MAX_OUTPUT).abs() / STEP_TEST_MAX_OUTPUT;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }
}
