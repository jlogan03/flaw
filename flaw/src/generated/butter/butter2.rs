//! Butterworth filter of order 2.
//! Region of validity: cutoff ratio from 1.00e-03 to 4.00e-01 .
//! This file is autogenerated.
#![allow(clippy::style)]

use crate::SisoIirFilter;

/// Minimum tabulated cutoff ratio
#[allow(dead_code)]
pub const MIN_CUTOFF_RATIO: f64 = 0.001;

/// Maximum tabulated cutoff ratio
#[allow(dead_code)]
pub const MAX_CUTOFF_RATIO: f64 = 0.4;

/// Initialise a Butterworth filter of order 2 by interpolating the coefficients from stored tables.
/// Cutoff ratio is the dimensionless ratio of the cutoff frequency to the sampling frequency.
/// Region of validity: cutoff ratio from 1.00e-03 to 4.00e-01
pub fn butter2(cutoff_ratio: f64) -> Result<SisoIirFilter<2>, &'static str> {
    let avals = &[&AVALS[0][..], &AVALS[1][..]];
    let cvals = &[&CVALS[0][..], &CVALS[1][..]];
    SisoIirFilter::new_interpolated(cutoff_ratio, &LOG10_CUTOFF_RATIOS, avals, cvals, &DVALS)
}

/// [dimensionless] Log base-10 of cutoff ratios, to improve float precision during interpolation
#[rustfmt::skip]
static LOG10_CUTOFF_RATIOS: [f64; 100] = [-3.0, -2.973716565744162, -2.947433131488324, -2.921149697232486, -2.894866262976648, -2.86858282872081, -2.842299394464972, -2.816015960209134, -2.789732525953296, -2.763449091697458, -2.73716565744162, -2.710882223185782, -2.684598788929944, -2.658315354674106, -2.632031920418268, -2.60574848616243, -2.579465051906592, -2.553181617650754, -2.526898183394916, -2.500614749139078, -2.4743313148832398, -2.448047880627402, -2.4217644463715637, -2.395481012115726, -2.3691975778598877, -2.34291414360405, -2.316630709348212, -2.290347275092374, -2.2640638408365357, -2.237780406580698, -2.21149697232486, -2.185213538069022, -2.1589301038131836, -2.132646669557346, -2.106363235301508, -2.08007980104567, -2.0537963667898316, -2.027512932533994, -2.001229498278156, -1.9749460640223178, -1.9486626297664797, -1.9223791955106417, -1.8960957612548037, -1.8698123269989657, -1.8435288927431277, -1.8172454584872897, -1.7909620242314517, -1.7646785899756137, -1.7383951557197759, -1.7121117214639379, -1.6858282872080999, -1.6595448529522618, -1.6332614186964238, -1.6069779844405858, -1.5806945501847478, -1.5544111159289098, -1.5281276816730718, -1.5018442474172338, -1.4755608131613958, -1.4492773789055577, -1.4229939446497197, -1.3967105103938817, -1.3704270761380437, -1.3441436418822057, -1.3178602076263677, -1.2915767733705297, -1.2652933391146917, -1.2390099048588537, -1.2127264706030156, -1.1864430363471776, -1.1601596020913396, -1.1338761678355016, -1.1075927335796636, -1.0813092993238256, -1.0550258650679876, -1.0287424308121496, -1.0024589965563115, -0.9761755623004738, -0.9498921280446355, -0.9236086937887977, -0.8973252595329595, -0.8710418252771217, -0.8447583910212835, -0.8184749567654457, -0.7921915225096074, -0.7659080882537697, -0.7396246539979314, -0.7133412197420936, -0.6870577854862554, -0.6607743512304176, -0.6344909169745794, -0.6082074827187416, -0.5819240484629034, -0.5556406142070656, -0.5293571799512273, -0.5030737456953895, -0.47679031143955175, -0.4505068771837135, -0.4242234429278757, -0.3979400086720376];

/// State-Space `D` 1x1 matrix
#[rustfmt::skip]
static DVALS: [f64; 100] = [9.825916820482034e-06, 1.1087150412788638e-05, 1.2510058263239782e-05, 1.4115322325479205e-05, 1.592626180963786e-05, 1.796916779788461e-05, 2.0273679984977082e-05, 2.2873210781532734e-05, 2.5805422656063134e-05, 2.9112765305368447e-05, 3.2843080040043865e-05, 3.705027966178098e-05, 4.17951131019404e-05, 4.7146025197536446e-05, 5.3180123213350725e-05, 5.9984263090454736e-05, 6.765626992618364e-05, 7.630630888390017e-05, 8.605842460926807e-05, 9.705226931091074e-05, 0.00010944504196621286, 0.00012341366365762777, 0.0001391572168523801, 0.00015689967953093638, 0.00017689298847002, 0.00019942046971769195, 0.00022480067838491743, 0.0002533916943436551, 0.0002855959252889097, 0.00032186547391311086, 0.00036270813167537714, 0.0004086940678425264, 0.0004604632891458456, 0.0005187339525450508, 0.0005843116212189892, 0.000658099562003302, 0.0007411101910495238, 0.0008344777834560534, 0.0009394725719716005, 0.0010575163695304129, 0.0011901998602590067, 0.0013393017135851497, 0.0015068096860431982, 0.001694943885137762, 0.0019061823789992753, 0.002143289344306371, 0.002409345952792225, 0.002707784203293268, 0.003042423911407351, 0.003417513072049895, 0.0038377718111634565, 0.004308440141185702, 0.004835329730275838, 0.005424879887461927, 0.006084217954622211, 0.006821224281527818, 0.007644601942250824, 0.00856395133058964, 0.009589849749710331, 0.010733936088423553, 0.012009000655594912, 0.013429080228196985, 0.015009558361642354, 0.016767271018845058, 0.018720617604160734, 0.02088967754917239, 0.023296332700823073, 0.02596439592308522, 0.02891974655898164, 0.032190473732074766, 0.03580702892195311, 0.0398023898588235, 0.0442122385869304, 0.049075157592453324, 0.05443284923666192, 0.0603303854510001, 0.06681649682704527, 0.07394391298489253, 0.0817697695757693, 0.09035610166235093, 0.09977044878241761, 0.11008660408991189, 0.12138554906624446, 0.1337566270779312, 0.14729902447272972, 0.162123648301674, 0.1783555170558988, 0.19613681779469241, 0.21563083375484057, 0.23702701688512964, 0.2605475784610542, 0.286456110928736, 0.31506895465004997, 0.3467703129548376, 0.38203254036297607, 0.42144364393295314, 0.46574493280111773, 0.5158830344086691, 0.5730822637195863, 0.6389455251590224];

/// State-Space `A` matrix, first row
#[rustfmt::skip]
static AVALS: [[f64; 100]; 2] = [[1.9911142922016536, 1.9905599325376073, 1.9899709886749763, 1.9893453032619322, 1.9886805844137065, 1.9879743973313582, 1.9872241554000196, 1.9864271107346145, 1.985580344139128, 1.984680754443502, 1.983725047180113, 1.9827097225595753, 1.9816310627032694, 1.9804851180875662, 1.9792676931521496, 1.977974331022193, 1.9766002972913506, 1.9751405628096552, 1.9735897854174407, 1.9719422905633148, 1.9701920507410806, 1.968332663677261, 1.966357329197634, 1.964258824697895, 1.9620294791402622, 1.959661145494636, 1.9571451715397308, 1.954472368936639, 1.9516329804844357, 1.9486166454649645, 1.9454123629817646, 1.9420084531964703, 1.938392516364942, 1.9345513895751076, 1.9304711010890923, 1.9261368221939532, 1.921532816468367, 1.916642386377247, 1.911447817112686, 1.905930317608211, 1.90006995866434, 1.8938456081372677, 1.8872348631594713, 1.8802139793815622, 1.8727577972491314, 1.864839665357043, 1.856431360956923, 1.8475030077317316, 1.8380229909944241, 1.8279578705158543, 1.8172722912399364, 1.8059288922012817, 1.7938882140210821, 1.7811086054196317, 1.7675461292465509, 1.7531544685897682, 1.737884833577955, 1.721685869533551, 1.7045035671585782, 1.6862811754353346, 1.6669591178890624, 1.6464749127780172, 1.624763097633663, 1.6017551583530407, 1.5773794627265738, 1.5515611978443595, 1.5242223102354868, 1.4952814468273627, 1.464653893830751, 1.432251509422256, 1.3979826445657086, 1.3617520444386908, 1.3234607206557023, 1.2830057817436242, 1.2402802060573377, 1.19517253744067, 1.1475664793403932, 1.097340357645242, 1.0443664160909862, 0.9885099004428943, 0.9296278785700803, 0.8675677326055805, 0.8021652461662265, 0.7332421934574906, 0.6606033171832251, 0.5840325574531908, 0.5032883629979747, 0.41809787736979737, 0.3281497446919729, 0.23308522038462734, 0.13248720166020841, 0.025866712972501527, -0.08735369718578544, -0.20786024337869474, -0.3364724310391757, -0.47417211015099114, -0.6221365148726473, -0.7817726558068178, -0.9547454589423299, -1.142980502539901], [-0.9911535958689355, -0.9906042811392585, -0.9900210289080292, -0.9894017645512341, -0.9887442894609451, -0.9880462740025496, -0.9873052501199595, -0.9865186035777406, -0.9856835658297523, -0.9847972055047233, -0.9838564195002732, -0.9828579236782224, -0.9817982431556772, -0.9806737021883564, -0.9794804136450029, -0.9782142680745549, -0.9768709223710553, -0.9754457880451908, -0.9739340191158777, -0.9723304996405585, -0.9706298309089455, -0.9688263183318915, -0.9669139580650434, -0.9648864234160187, -0.9627370510941422, -0.9604588273735067, -0.9580443742532704, -0.9554859357140136, -0.9527753641855913, -0.9499041073606169, -0.9468631955084661, -0.9436432294678404, -0.9402343695215254, -0.9366263253852879, -0.9328083475739682, -0.9287692204419664, -0.924497257232565, -0.9199802975110712, -0.9152057074005725, -0.9101603830863326, -0.9048307581053759, -0.8992028149916083, -0.8932621019036441, -0.8869937549221132, -0.8803825267651285, -0.8734128227342685, -0.8660687447680919, -0.8583341445449046, -0.8501926866400535, -0.8416279228040539, -0.8326233784845902, -0.8231626527660246, -0.8132295329421855, -0.8028081249694794, -0.7918830010650397, -0.7804393657158795, -0.7684632413469584, -0.7559416748559096, -0.7428629661574195, -0.7292169197890288, -0.7149951205114421, -0.7001912336908052, -0.6848013310802322, -0.668824242428421, -0.6522619331432168, -0.6351199080410491, -0.6174076410387791, -0.5991390305197036, -0.5803328800666775, -0.561013404350555, -0.5412107602535211, -0.5209616038739848, -0.5003096750034238, -0.4793064121134375, -0.4580116030039854, -0.43649407924467043, -0.41483246664857426, -0.3931160095848123, -0.3714454943940635, -0.34993430709229795, -0.3287096736997508, -0.307914148965228, -0.2877074424312045, -0.2682687017692154, -0.24979941507414416, -0.23252715065988672, -0.21671043122156997, -0.20264514854856683, -0.1906730797113353, -0.18119328792514583, -0.1746775155044252, -0.1716911566874457, -0.17292212141441457, -0.17922100844065597, -0.19165773041272846, -0.21160246558082108, -0.2408432163318235, -0.281759481827858, -0.33758359593601506, -0.41280159809618855]];

/// State-Space `C` vector
#[rustfmt::skip]
static CVALS: [[f64; 100]; 2] = [[3.921635705621048e-05, 4.424393820329213e-05, 4.991476953696039e-05, 5.6310894823178765e-05, 6.352477126239203e-05, 7.166058111931492e-05, 8.083570655495053e-05, 9.118238756904937e-05, 0.00010284958531020776, 0.00011600507562093215, 0.000130837800582064, 0.00014756050903252642, 0.00016641272059588229, 0.00018766405167577524, 0.00021161794622073722, 0.0002386158588391109, 0.0002690419431020857, 0.0003033283036426743, 0.00034196087697650453, 0.0003854860128711465, 0.0004345178356072823, 0.0004897464726366436, 0.0005519472469731823, 0.0006219909391727862, 0.0007008552349714377, 0.0007896374855574345, 0.000889568919029732, 0.0010020304538000226, 0.0011285702774636245, 0.0012709233678937788, 0.0014310331468660507, 0.0016110754702064907, 0.0018134851720327848, 0.002040985393805918, 0.0022966199412317533, 0.0025837889230508762, 0.0029062879348198487, 0.0032683510571741127, 0.003674697940875558, 0.004130585249115808, 0.004641862719802647, 0.0052150340954142515, 0.0058573231437334975, 0.006576744957378842, 0.007382182671248343, 0.00828346967221235, 0.009291477291742574, 0.01041820786645932, 0.011676892920332596, 0.013082096061744212, 0.014649819994855968, 0.016397616813658433, 0.01834470047459924, 0.020512060025850175, 0.02292257180442928, 0.025601108393469164, 0.028574641658679935, 0.03187233665442851, 0.03552563260631672, 0.03956830654108173, 0.04403651445076949, 0.0489688041538036, 0.05440609326105998, 0.06039160488362863, 0.06697075294668217, 0.07419096821912112, 0.08210145545090922, 0.09075287134803979, 0.1001969125241742, 0.1104858020559311, 0.12167166283025904, 0.1338057654814457, 0.14693763831592155, 0.16111402611600378, 0.17637768394085887, 0.19276599076624568, 0.21030936567976144, 0.2290294658903157, 0.24893714034797054, 0.2700301043833604, 0.29229028821041897, 0.31568079368035684, 0.3401423669802355, 0.36558925678396026, 0.39190427312999776, 0.4189327855446163, 0.4464752903225183, 0.4742780227834124, 0.5020208705540493, 0.5293015282380406, 0.555614376491757, 0.5803218998580847, 0.6026154712429576, 0.6214608642623766, 0.6355216631339495, 0.6430504659124949, 0.6417329362897528, 0.628462818921993, 0.5990168385525069, 0.5475887728761645], [8.692403115220801e-08, 1.041717482453162e-07, 1.2483750976774014e-07, 1.4959750944065066e-07, 1.792613928984892e-07, 2.1479850825811988e-07, 2.573692965572671e-07, 3.083628219957423e-07, 3.694416346909479e-07, 4.425953881267403e-07, 5.302049064854187e-07, 6.351187217054502e-07, 7.607444859624883e-07, 9.111581236028466e-07, 1.091234130645732e-06, 1.3068010754340185e-06, 1.5648271192075303e-06, 1.8736412818244234e-06, 2.2431972507828583e-06, 2.6853878005828573e-06, 3.214419388725232e-06, 3.847258264357903e-06, 4.6041615123297e-06, 5.5093089132116654e-06, 6.591554391162851e-06, 7.88531921836364e-06, 9.431653129928516e-06, 1.1279494171548461e-05, 1.3487163561847837e-05, 1.6124138225475347e-05, 1.9273151080324076e-05, 2.30326777992561e-05, 2.7519878787993592e-05, 3.287407672019358e-05, 3.926086336143741e-05, 4.687694482829563e-05, 5.595585211713673e-05, 6.67746639657741e-05, 7.966191215689664e-05, 9.500686551854455e-05, 0.0001132704184039371, 0.00013499784260629852, 0.0001608336987194809, 0.0001915392440771436, 0.00022801271970072944, 0.00027131294815946404, 0.00032268672774538043, 0.00038360056554733466, 0.00045577735226999495, 0.0005412386440648413, 0.0006423532798996144, 0.0007618931252836537, 0.0009030967921021546, 0.0010697422368239774, 0.001266229181582177, 0.0014976723298464924, 0.0017700063549015016, 0.002090103618359212, 0.0024659055196365274, 0.0029065682768110324, 0.0034226237846258414, 0.004026155975882939, 0.00473099281666324, 0.005552913682073994, 0.006509871376035918, 0.007622227465134853, 0.00891299888315333, 0.010408112921698196, 0.012136666747609437, 0.014131186475986386, 0.016427879576683048, 0.01906687299995318, 0.022092427868329417, 0.025553119882913, 0.02950197270170413, 0.033996529403089754, 0.039098844635465445, 0.04487537697918499, 0.051396757089209016, 0.05873740183557493, 0.06697493711827142, 0.07618938107909468, 0.08646202319628779, 0.09787391037870552, 0.11050381431844979, 0.12442549830750016, 0.1397040160439689, 0.1563906432168438, 0.17451583860208217, 0.19407931236862394, 0.21503577478478297, 0.23727412990319408, 0.26058656262014135, 0.28462178776978975, 0.308813050733199, 0.33226512977337463, 0.35357342519504753, 0.3705280979498995, 0.37961909236597663, 0.3751877912769695]];

#[cfg(feature = "std")]
#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    static CUTOFF_TEST_INPUT: [f32; 101] = [0.0, 0.5877852522924732, -0.9510565162951536, 0.9510565162951535, -0.5877852522924728, -4.898587196589413e-16, 0.5877852522924736, -0.9510565162951538, 0.9510565162951533, -0.5877852522924725, -9.797174393178826e-16, 0.587785252292474, -0.951056516295154, 0.9510565162951532, -0.587785252292472, -1.4695761589768238e-15, 0.5877852522924744, -0.9510565162951541, 0.951056516295153, -0.5877852522924716, -1.959434878635765e-15, 0.5877852522924748, -0.9510565162951542, 0.9510565162951529, -0.5877852522924712, -2.4492935982947065e-15, 0.5877852522924751, -0.9510565162951544, 0.9510565162951528, -0.5877852522924708, -2.9391523179536475e-15, 0.5877852522924756, -0.9510565162951545, 0.9510565162951525, -0.5877852522924705, -3.429011037612589e-15, 0.587785252292476, -0.9510565162951546, 0.9510565162951524, -0.58778525229247, -3.91886975727153e-15, 0.5877852522924764, -0.9510565162951549, 0.9510565162951523, -0.5877852522924697, -4.408728476930472e-15, 0.5877852522924768, -0.951056516295155, 0.9510565162951521, -0.5877852522924693, -4.898587196589413e-15, 0.5877852522924887, -0.9510565162951552, 0.9510565162951563, -0.5877852522924688, -1.9599300631450357e-14, 0.5877852522924776, -0.9510565162951509, 0.9510565162951519, -0.5877852522924569, -5.878304635907295e-15, 0.5877852522924665, -0.9510565162951554, 0.9510565162951473, -0.587785252292468, 7.842691359635767e-15, 0.5877852522924784, -0.95105651629516, 0.9510565162951515, -0.5877852522924791, -6.858022075225178e-15, 0.5877852522924902, -0.9510565162951558, 0.9510565162951558, -0.5877852522924673, -2.1558735510086122e-14, 0.5877852522924791, -0.9510565162951515, 0.9510565162951512, -0.5877852522924554, -7.83773951454306e-15, 0.587785252292468, -0.9510565162951561, 0.9510565162951466, -0.5877852522924665, 5.883256481000002e-15, 0.5877852522924799, -0.9510565162951606, 0.9510565162951509, -0.5877852522924776, -8.817456953860943e-15, 0.5877852522924918, -0.9510565162951563, 0.9510565162951552, -0.5877852522924657, -2.3518170388721888e-14, 0.5877852522924807, -0.9510565162951521, 0.9510565162951506, -0.5877852522924538, -9.797174393178826e-15];
    static CUTOFF_TEST_OUTPUT: [f32; 101] = [0.0, 0.3755627567067428, -0.28580870024257504, -0.06046968287379373, 0.41920850528282144, -0.5976373744638563, 0.510037925685062, -0.19280553319986227, -0.2222820541365811, 0.5657650347226963, -0.6983522246845361, 0.5646542662393135, -0.21365569445033644, -0.22099643901543709, 0.572902581590217, -0.7070410055666215, 0.5716389826249801, -0.2180523464612913, -0.2188544535767145, 0.572269278973428, -0.7072013680356273, 0.572083702132682, -0.2184944543042131, -0.21853271385571882, 0.5720840395695259, -0.7071224576796634, 0.5720699762563272, -0.21851134021620805, -0.21850774752384672, 0.5720624740704332, -0.7071081148763677, 0.5720624849842979, -0.21850869857046068, -0.21850767446436514, 0.572061300089284, -0.7071068031978748, 0.572061470382649, -0.21850808036353608, -0.21850796223364433, 0.5720613738071529, -0.7071067686643433, 0.5720614004806417, -0.21850801472240167, -0.2185080084045209, 0.5720613994827992, -0.7071067789516949, 0.5720614016399359, -0.21850801180081691, -0.2185080122223937, 0.5720614026405185, -0.7071067809848823, 0.5720614026603253, -0.2185080121277843, -0.2185080122698847, 0.5720614028297779, -0.7071067811816042, 0.5720614028070222, -0.2185080122142699, -0.21850801223159466, 0.5720614028217196, -0.707106781188178, 0.5720614028178809, -0.2185080122239757, -0.2185080122249995, 0.572061402818172, -0.7071067811868542, 0.5720614028178501, -0.21850801222446586, -0.21850801222442184, 0.5720614028177095, -0.707106781186585, 0.5720614028177153, -0.2185080122244143, -0.2185080122244094, 0.5720614028176898, -0.7071067811865591, 0.5720614028176751, -0.2185080122244002, -0.21850801222441374, 0.5720614028176929, -0.7071067811865401, 0.57206140281767, -0.21850801222441074, -0.2185080122244154, 0.5720614028176836, -0.7071067811865369, 0.5720614028176885, -0.2185080122244119, -0.21850801222441707, 0.5720614028176816, -0.7071067811865546, 0.5720614028176911, -0.21850801222439942, -0.21850801222441696, 0.5720614028176918, -0.7071067811865581, 0.5720614028176728, -0.21850801222439747, -0.21850801222441585, 0.5720614028176947, -0.7071067811865402];
    const STEP_TEST_MIN_OUTPUT: f32 = 1.000000000001043;
    const STEP_TEST_MAX_OUTPUT: f32 = 1.1865342980351867;

    #[test]
    fn test() {
        let order = 2;
        println!("order {order}");
        let mut filter = butter2(0.4).unwrap();
        let out = (0..CUTOFF_TEST_INPUT.len()).map(|i| {filter.update(CUTOFF_TEST_INPUT[i])}).collect::<Vec<f32>>();
        // Check overall match to reference output to catch phase error, etc
        (0..CUTOFF_TEST_INPUT.len()).for_each(|i| { let expected = CUTOFF_TEST_OUTPUT[i]; let rel_err = (out[i] - expected).abs() / expected.abs().max(1e-4); assert!(rel_err < 0.05); });
        // Check approximate attenuation at cutoff frequency; should be -3dB or 1/sqrt(2) magnitude
        let maxmag = out.iter().fold(0.0_f32, |a, b| a.abs().max(b.abs()));
        let attenuation_rel_err = (maxmag - 0.707).abs() / 0.707;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);
        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter2(MIN_CUTOFF_RATIO).unwrap();
        (0..9999).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        let step_min_rel_err = (step_min_final - STEP_TEST_MIN_OUTPUT).abs() / STEP_TEST_MIN_OUTPUT;
        println!("order {order} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);
        let mut filtermax = butter2(MAX_CUTOFF_RATIO).unwrap();
        (0..1).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - STEP_TEST_MAX_OUTPUT).abs() / STEP_TEST_MAX_OUTPUT;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }
}
