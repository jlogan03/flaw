//! Butterworth filter of order 5.
//! Region of validity: cutoff ratio from 5.62e-02 to 4.00e-01 .
//! This file is autogenerated.
#![allow(clippy::style)]

use crate::SisoIirFilter;

/// Minimum tabulated cutoff ratio
#[allow(dead_code)]
pub const MIN_CUTOFF_RATIO: f64 = 0.05623413251903491;

/// Maximum tabulated cutoff ratio
#[allow(dead_code)]
pub const MAX_CUTOFF_RATIO: f64 = 0.4;

/// Initialise a Butterworth filter of order 5 by interpolating the coefficients from stored tables.
/// Cutoff ratio is the dimensionless ratio of the cutoff frequency to the sampling frequency.
/// Region of validity: cutoff ratio from 5.62e-02 to 4.00e-01
pub fn butter5(cutoff_ratio: f64) -> Result<SisoIirFilter<5>, &'static str> {
    let avals = &[
        &AVALS[0][..],
        &AVALS[1][..],
        &AVALS[2][..],
        &AVALS[3][..],
        &AVALS[4][..],
    ];
    let cvals = &[
        &CVALS[0][..],
        &CVALS[1][..],
        &CVALS[2][..],
        &CVALS[3][..],
        &CVALS[4][..],
    ];
    SisoIirFilter::new_interpolated(cutoff_ratio, &LOG10_CUTOFF_RATIOS, avals, cvals, &DVALS)
}

/// [dimensionless] Log base-10 of cutoff ratios, to improve float precision during interpolation
#[rustfmt::skip]
static LOG10_CUTOFF_RATIOS: [f64; 100] = [-1.25, -1.2413933334209297, -1.2327866668418594, -1.224180000262789, -1.2155733336837187, -1.2069666671046484, -1.1983600005255781, -1.1897533339465076, -1.1811466673674373, -1.172540000788367, -1.1639333342092968, -1.1553266676302263, -1.146720001051156, -1.1381133344720857, -1.1295066678930155, -1.1209000013139452, -1.1122933347348747, -1.1036866681558044, -1.0950800015767341, -1.0864733349976639, -1.0778666684185936, -1.069260001839523, -1.0606533352604528, -1.0520466686813825, -1.043440002102312, -1.0348333355232417, -1.0262266689441715, -1.0176200023651012, -1.009013335786031, -1.0004066692069604, -0.9918000026278901, -0.9831933360488199, -0.9745866694697496, -0.9659800028906792, -0.9573733363116088, -0.9487666697325385, -0.9401600031534683, -0.9315533365743979, -0.9229466699953276, -0.9143400034162572, -0.9057333368371869, -0.8971266702581167, -0.8885200036790463, -0.8799133370999759, -0.8713066705209056, -0.8627000039418353, -0.8540933373627649, -0.8454866707836945, -0.8368800042046243, -0.828273337625554, -0.8196666710464837, -0.8110600044674133, -0.8024533378883429, -0.7938466713092727, -0.7852400047302024, -0.776633338151132, -0.7680266715720616, -0.7594200049929913, -0.7508133384139211, -0.7422066718348507, -0.7336000052557804, -0.72499333867671, -0.7163866720976397, -0.7077800055185693, -0.6991733389394991, -0.6905666723604288, -0.6819600057813584, -0.6733533392022881, -0.6647466726232177, -0.6561400060441475, -0.6475333394650771, -0.6389266728860068, -0.6303200063069364, -0.6217133397278661, -0.6131066731487957, -0.6045000065697255, -0.5958933399906552, -0.5872866734115848, -0.5786800068325145, -0.5700733402534441, -0.5614666736743739, -0.5528600070953036, -0.5442533405162332, -0.5356466739371628, -0.5270400073580925, -0.5184333407790221, -0.5098266741999519, -0.5012200076208816, -0.4926133410418112, -0.4840066744627409, -0.47540000788367054, -0.46679334130460026, -0.4581866747255299, -0.4495800081464596, -0.4409733415673892, -0.43236667498831893, -0.42376000840924855, -0.41515334183017827, -0.406546675251108, -0.3979400086720376];

/// State-Space `D` 1x1 matrix
#[rustfmt::skip]
static DVALS: [f64; 100] = [0.00010200173380017545, 0.00011155755977679115, 0.00012198962436515001, 0.0001333762277085027, 0.00014580245275948816, 0.00015936073142622627, 0.0001741514559613601, 0.00019028363903946056, 0.00020787562621904922, 0.0002270558647543096, 0.0002479637330096711, 0.00027075043503942256, 0.00029557996522590297, 0.00032263014822549716, 0.00035209375985353154, 0.0003841797349494436, 0.0004191144687046418, 0.0004571432184098989, 0.0004985316130899056, 0.0005435672790428272, 0.0005925615898961423, 0.0006458515504303865, 0.0007038018241143105, 0.0007668069150431365, 0.0008352935157815311, 0.0009097230334907522, 0.0009905943076718806, 0.0010784465338916733, 0.0011738624089829939, 0.00127747151443729, 0.0013899539560429987, 0.001512044279282924, 0.0016445356815989133, 0.0017882845443787391, 0.0019442153094347233, 0.0021133257268455205, 0.0022966925033428886, 0.002495477382968509, 0.00271093369452911, 0.0029444134034721783, 0.003197374709223804, 0.0034713902328145163, 0.0037681558438119266, 0.004089500180231691, 0.004437394920266975, 0.0048139658704263125, 0.00522150494107343, 0.005662483087503895, 0.006139564302666913, 0.0066556207565540235, 0.007213749187252935, 0.007817288659844874, 0.008469839821867263, 0.009175285798154452, 0.009937814883717413, 0.010761945211169148, 0.011652551589324298, 0.012614894732315515, 0.013654653124243896, 0.014777957793437541, 0.015991430303321293, 0.017302224304266478, 0.018718071033246794, 0.020247329196438933, 0.021899039724955098, 0.023682985956716227, 0.025609759869267088, 0.02769083507051088, 0.029938647348549015, 0.032366683689987664, 0.03498958080049848, 0.037823234304779886, 0.04088491996853519, 0.044193428476419515, 0.04776921552157453, 0.05163456921966562, 0.05581379715959837, 0.06033343575186427, 0.06522248494282779, 0.07051267184008685, 0.07623874735343536, 0.08243882061377286, 0.08915473670771693, 0.09643250418201478, 0.10432277985751705, 0.11288141978224847, 0.12217010668958353, 0.13225706616324565, 0.14321788591060036, 0.15513645518998082, 0.16810604462669698, 0.1822305505111344, 0.19762593235818335, 0.21442187821835176, 0.2327637392185457, 0.2528147833944003, 0.27475882946543906, 0.29880333432557044, 0.3251830243552488, 0.35416418109342973];

/// State-Space `A` matrix, first row
#[rustfmt::skip]
static AVALS: [[f64; 100]; 5] = [[3.8582729539197853, 3.8354872729103455, 3.812249377530346, 3.7885504311741762, 3.7643814314722768, 3.7397332074610943, 3.7145964167088543, 3.6889615423955533, 3.6628188903453687, 3.6361585860094783, 3.6089705713970446, 3.5812446019518624, 3.5529702433718637, 3.5241368683683914, 3.4947336533617617, 3.4647495751093036, 3.4341734072615977, 3.402993716842219, 3.371198860645774, 3.3387769815484774, 3.3057160047249483, 3.272003633764236, 3.2376273466774195, 3.202574391788366, 3.1668317834984405, 3.1303862979150656, 3.0932244683331214, 3.0553325805571583, 3.016696668051325, 2.9773025069027628, 2.9371356105829896, 2.8961812244904808, 2.8544243202562787, 2.811849589792971, 2.768441439065806, 2.7241839815630926, 2.679061031441228, 2.63305609631794, 2.5861523696853372, 2.538332722912374, 2.489579696804264, 2.4398754926841324, 2.38920196296001, 2.3375406011378796, 2.2848725312391243, 2.2311784965782793, 2.1764388478544845, 2.120633530507559, 2.0637420712870727, 2.0057435639803103, 1.9466166542426198, 1.8863395234711846, 1.824889871661128, 1.7622448991806836, 1.698381287400379, 1.6332751781095742, 1.5669021516525046, 1.4992372037151913, 1.4302547206943659, 1.359928453580022, 1.2882314902844163, 1.215136226352546, 1.140614333992505, 1.0646367293687702, 0.9871735381078373, 0.9081940589737392, 0.8276667256814833, 0.7455590668294061, 0.6618376639475287, 0.576468107678643, 0.489414952132589, 0.4006416674828042, 0.3101105909083524, 0.217782876025308, 0.12361844099939956, 0.027575915588579306, -0.07038741256937464, -0.17031565203114807, -0.27225436150037946, -0.37625060306921904, -0.4823529940379375, -0.5906117564611613, -0.7010787634080646, -0.8138075807399595, -0.9288535029997015, -1.0462735817717204, -1.1661266446074263, -1.2884733023180204, -1.413375942115319, -1.540898703733797, -1.6711074352981274, -1.8040696253184094, -1.9398543068137841, -2.078531929203735, -2.2201741932950094, -2.3648538444713245, -2.5126444191198427, -2.6636199394784934, -2.8178545525615943, -2.9754221097456828], [-6.058101770819958, -5.991614961127277, -5.924278399614318, -5.856093226102329, -5.7870612573166875, -5.7171850222028935, -5.646467798579294, -5.574913651168228, -5.502527471048573, -5.429315016574163, -5.3552829558040935, -5.2804389104925, -5.204791501687283, -5.1283503969891875, -5.051126359524709, -4.973131298688727, -4.894378322715191, -4.814881793137093, -4.7346573811999235, -4.653722126296277, -4.572094496492962, -4.489794451226059, -4.406843506244005, -4.323264800883709, -4.239083167770396, -4.1543252050379404, -4.06901935117333, -3.9831959625964757, -3.8968873940949664, -3.8101280822425925, -3.7229546319408606, -3.635405906233942, -3.5475231195602595, -3.4593499346177254, -3.370932563035092, -3.2823198700589473, -3.1935634834845965, -3.104717907079972, -3.015840638774554, -2.926992293910546, -2.8382367338815775, -2.7496412005147066, -2.661276456585549, -2.573216932893404, -2.4855408823643437, -2.398330541695217, -2.3116723011010096, -2.2256568827824372, -2.1403795287902856, -2.0559401990285657, -1.9724437802104777, -1.8900003066599051, -1.808725193937656, -1.7287394863662513, -1.6501701196307996, -1.5731501997469703, -1.497819299811338, -1.4243237760852538, -1.3528171051120341, -1.28346024372977, -1.21642201401955, -1.1518795154228059, -1.0900185664732063, -1.031034178819419, -0.975131066466949, -0.9225241934415873, -0.8734393633757961, -0.8281138548442442, -0.7867971066277196, -0.749751457467728, -0.7172529452890064, -0.6895921713160001, -0.6670752349935324, -0.6500247461431073, -0.6387809213454225, -0.6337027721374145, -0.63516939324818, -0.6435813597711898, -0.6593622428776593, -0.6829602544132134, -0.7148500314801429, -0.7555345728804821, -0.8055473400667206, -0.865454535997929, -0.9358575760036316, -1.0173957653821022, -1.1107491989596605, -1.2166418981567109, -1.3358452011724715, -1.4691814216233998, -1.6175277902373093, -1.7818206928769367, -1.9630602160742376, -2.1623150081953595, -2.3807274600854638, -2.619519203281491, -2.879996916313239, -3.1635584198930937, -3.4716990295662042, -3.806018119320411], [4.823280869710185, 4.748782768529209, 4.673774851507496, 4.598275835908028, 4.522305420670898, 4.445884303384452, 4.369034196334396, 4.291777841524347, 4.214139024555109, 4.136142587243691, 4.057814438856416, 3.979181565823619, 3.900272039796239, 3.8211150238972147, 3.7417407770127573, 3.66218065596056, 3.582467115363446, 3.5026337050482366, 3.422715064780283, 3.342746916134455, 3.262766051293237, 3.1828103185517578, 3.102918604298439, 3.023130811227833, 2.943487832529735, 2.8640315217850403, 2.784804658284545, 2.705850907471538, 2.6272147761925484, 2.5489415624228675, 2.4710772991143988, 2.393668691792536, 2.3167630495063793, 2.2404082087119033, 2.1646524496408768, 2.0895444046789025, 2.0151329582434143, 1.9414671376168566, 1.8685959941506542, 1.7965684742117793, 1.725433279195189, 1.6552387138711215, 1.5860325222760967, 1.5178617102889347, 1.4507723539579263, 1.3848093925608835, 1.3200164052852623, 1.256435370309499, 1.194106404947429, 1.133067485383695, 1.0733541443772279, 1.014999145139973, 0.9580321294066931, 0.90247923749572, 0.8483626979169979, 0.79570038380866, 0.7445053331726491, 0.6947852295285964, 0.6465418392078286, 0.5997704010597158, 0.5544589638333094, 0.5105876659202018, 0.4681279514903535, 0.42704171631033716, 0.38728037569084073, 0.3487838460529058, 0.31147943051407784, 0.2752805976573982, 0.24008564123640014, 0.2057762069629332, 0.1722156706928279, 0.13924735023391216, 0.10669253061308634, 0.07434827990966418, 0.041985029639399575, 0.009343890098084095, -0.023866333023784055, -0.057972459142681496, -0.09334063418056528, -0.1303804628347066, -0.16954962881219998, -0.21135906830921444, -0.25637877143473753, -0.3052442971217608, -0.3586640995544595, -0.4174277785194274, -0.48241538264957196, -0.5546079135981075, -0.6350992011270334, -0.7251093443353755, -0.8259999432494649, -0.939291378258057, -1.0666824329468811, -1.2100725993459096, -1.37158745403327, -1.5536075495016477, -1.7588013281692654, -1.9901626367350485, -2.251053496291228, -2.5452528683304663], [-1.943091498433082, -1.9053004841751737, -1.8674468106727737, -1.8295449886762822, -1.791609921237361, -1.753656896859613, -1.7157015815836791, -1.677760009968039, -1.639848574927902, -1.6019840163959969, -1.5641834087707116, -1.5264641471189704, -1.4888439321034923, -1.4513407536066643, -1.4139728730261565, -1.3767588042207464, -1.3397172930884433, -1.3028672957631744, -1.2662279554207638, -1.2298185776899904, -1.1936586046699975, -1.1577675875613176, -1.1221651579243945, -1.086870997586599, -1.051904807226572, -1.0172862736731196, -0.9830350359650764, -0.9491706502284563, -0.9157125534379243, -0.8826800261412051, -0.8500921542375944, -0.8179677899152293, -0.7863255118664728, -0.7551835849165427, -0.7245599192177025, -0.6944720291799144, -0.664936992329018, -0.6359714083054833, -0.6075913582406736, -0.5798123647736791, -0.5526493530003997, -0.5261166126778707, -0.5002277620414268, -0.47499571363032533, -0.4504326425596764, -0.4265499577233184, -0.40335827646442424, -0.38086740330887536, -0.35908631342168373, -0.33802314152011226, -0.3176851770598297, -0.2980788666038957, -0.2792098243903792, -0.2610828522348535, -0.24370197004133587, -0.2270704583520829, -0.21119091454630312, -0.19606532450401531, -0.1816951517884291, -0.16808144667355235, -0.15522497765934407, -0.1431263884818997, -0.13178638404940082, -0.12120594922594723, -0.1113866049569225, -0.10233070689532159, -0.09404179246537923, -0.08652498320783972, -0.07978745031417245, -0.07383895250341696, -0.068692456859137, -0.06436485496595945, -0.06087778871424755, -0.05825860253643482, -0.05654144167025728, -0.05576851939809917, -0.055991580191355576, -0.05727359041951893, -0.05969069391794378, -0.06333447643137328, -0.06831459098856021, -0.07476180589331041, -0.08283154857800575, -0.0927080324734964, -0.10460907081671861, -0.11879170157693858, -0.13555877221520418, -0.15526666276889733, -0.1783343619785846, -0.20525415534743904, -0.23660423801595087, -0.2730636314872378, -0.31542986448478827, -0.36463997824883543, -0.42179554001456265, -0.48819250110477497, -0.5653569253797899, -0.6550878490230853, -0.759508825558258, -0.8811300754378363], [0.3163753901414638, 0.30907556195003705, 0.30179731326956255, 0.2945439084097348, 0.28731864792257, 0.28012486481132026, 0.2729659205289593, 0.26584520076710466, 0.2587661110369863, 0.25173207204485326, 0.24474651486503565, 0.23781287591472747, 0.23093459173544367, 0.22411509358702977, 0.21735780186103326, 0.21066612032122678, 0.2040434301800427, 0.19749308402069488, 0.19101839957575395, 0.18462265337396516, 0.17830907426809833, 0.17208083685761064, 0.16594105482088295, 0.15989277417272782, 0.15393896646378305, 0.14808252193925023, 0.14232624267523983, 0.13667283571170283, 0.13112490620156125, 0.125684950596174, 0.12035534988769081, 0.11513836292910132, 0.110036119852909, 0.10505061560927478, 0.10018370364420066, 0.09543708973781012, 0.09081232602199948, 0.08631080519566649, 0.08193375495430426, 0.07768223264896207, 0.07355712018736248, 0.06955911918725885, 0.06568874638888736, 0.06194632932950098, 0.05833200227842622, 0.05484570242573074, 0.05148716631133696, 0.048255926474129685, 0.04515130829212667, 0.04217242697494329, 0.03931818465836597, 0.036587267537607146, 0.0339781429604618, 0.031489056383758936, 0.029118028075801536, 0.026862849423406342, 0.02472107867410996, 0.022690035911384927, 0.020766797022464173, 0.018948186373583312, 0.017230767854887037, 0.01561083389543068, 0.014084391975851231, 0.012647148080212814, 0.011294486426630728, 0.010021444695344718, 0.008822683829067611, 0.007692451308931241, 0.006624536604394412, 0.005612217249963032, 0.0046481937067749, 0.0037245108122866146, 0.002832463193215111, 0.0019624814991457107, 0.0011039956864959237, 0.00024527081955066145, -0.00062679007445353, -0.0015268826951173665, -0.0024715856939406832, -0.0034797021342663293, -0.004572669991090858, -0.005775056096563306, -0.007115151159413761, -0.008625687491325485, -0.01034470606603425, -0.012316605781762987, -0.014593415634809793, -0.017236340382125133, -0.020317642745803735, -0.02392294103937469, -0.028154021253452427, -0.03313228841566027, -0.039003015142176294, -0.04594058799341676, -0.05415500756515642, -0.06389997026156988, -0.07548295391191268, -0.08927785328853179, -0.10574087539067953, -0.1254306221553556]];

/// State-Space `C` vector
#[rustfmt::skip]
static CVALS: [[f64; 100]; 5] = [[0.0009035591997750197, 0.0009856653996047733, 0.0010750029913769539, 0.0011721837035359465, 0.0012778683096283717, 0.0013927702764110787, 0.0015176596540854986, 0.0016533672217609499, 0.0018007889016527527, 0.0019608904559017377, 0.0021347124802540123, 0.0023233757091581635, 0.0025280866471140383, 0.002750143541336119, 0.0029909427109664688, 0.0032519852481789078, 0.0035348841065472626, 0.0038413715919954107, 0.004173307271494118, 0.004532686314405265, 0.00492164828098565, 0.005342486372032423, 0.005797657152965496, 0.006289790764779068, 0.006821701633234764, 0.00739639968639094, 0.008017102089041573, 0.008687245500846536, 0.009410498862844669, 0.010190776714617461, 0.011032253041579588, 0.011939375648672066, 0.01291688105307974, 0.013969809884438163, 0.015103522776278851, 0.016323716727125348, 0.017636441903623578, 0.019048118851291318, 0.020565556070811834, 0.022195967909176117, 0.023946992705278038, 0.025826711118659786, 0.02784366455783432, 0.030006873610810707, 0.03232585636491291, 0.03481064648548849, 0.0374718109033835, 0.04032046693881221, 0.04336829866412055, 0.04662757227952211, 0.0501111502437004, 0.05383250386467284, 0.057805724014853975, 0.06204552958709491, 0.0665672732547417, 0.07138694403742352, 0.07652116610417554, 0.08198719316521579, 0.08780289771161338, 0.09398675425628805, 0.10055781560803342, 0.10753568106992406, 0.11494045529144517, 0.12279269631634421, 0.1311133511512235, 0.1399236769282293, 0.1492451454428208, 0.15909932850745145, 0.16950776116565763, 0.18049177934853872, 0.1920723280151077, 0.20426973518535924, 0.21710344653335908, 0.23059171433711095, 0.24475123355841402, 0.25959671662058226, 0.2751403970302558, 0.2913914503099632, 0.3083553187205612, 0.3260329238965796, 0.3444197487195454, 0.3635047664255773, 0.38326919097556683, 0.4036850179770124, 0.42471331977426385, 0.4463022515201922, 0.46838471686266225, 0.49087563202197815, 0.5136687151263428, 0.5366327132458063, 0.5596069620392528, 0.5823961515734822, 0.6047641457678052, 0.6264266708950796, 0.647042649144864, 0.6662039045225687, 0.6834229078669593, 0.6981181523356043, 0.7095966561810585, 0.7170329705617839], [0.0004020804538402056, 0.00044716565358243917, 0.0004971957470479774, 0.0005526986534781826, 0.000614256802008701, 0.000682512527424944, 0.0007581739714520813, 0.0008420215335195492, 0.0009349149183587489, 0.0010378008314312907, 0.0011517213770524625, 0.0012778232181792867, 0.001417367561182227, 0.001571741033522063, 0.0017424675271149796, 0.0019312210852954188, 0.002139839916682125, 0.0023703416249210734, 0.0026249397492214303, 0.0029060617168160043, 0.0032163683149641514, 0.0035587747968657685, 0.003936473742862271, 0.004352959805551148, 0.004812056474918066, 0.005317945007273317, 0.005875195669639776, 0.006488801459243257, 0.00716421446586216, 0.007907385052950407, 0.008724804041595181, 0.009623548089436853, 0.010611328464575197, 0.01169654342211291, 0.012888334398222375, 0.014196646243323386, 0.01563229172196021, 0.017207020512069696, 0.018933592940306967, 0.020825858692671798, 0.022898840740535113, 0.025168824720934027, 0.027653454006237334, 0.030371830691474096, 0.033344622717150316, 0.036594177330541294, 0.04014464106839279, 0.04402208641769676, 0.04825464527754927, 0.0528726493026519, 0.05790877715613391, 0.06339820863409293, 0.06937878554424488, 0.07589117912357943, 0.08297906366164143, 0.09068929585307478, 0.09907209923070576, 0.10818125282310553, 0.11807428293159034, 0.1288126566229818, 0.14046197517659356, 0.15309216529532962, 0.1667776653776646, 0.18159760353305246, 0.19763596328795344, 0.2149817320491542, 0.23372902633625126, 0.2539771865310118, 0.27583083237530415, 0.2993998686299114, 0.32479942812139956, 0.3521497367773719, 0.3815758820896495, 0.413207462637617, 0.4471780917129259, 0.4836247225440328, 0.5226867559192445, 0.5645048828977902, 0.6092196054703205, 0.6569693661016074, 0.7078882025887437, 0.762102827016531, 0.819729006067919, 0.8808670936581838, 0.9455965346957542, 1.0139691193457057, 1.086000718753564, 1.1616611736109674, 1.2408619334902615, 1.323440954118177, 1.4091442473764142, 1.497603339336249, 1.5883077181049055, 1.6805711368265377, 1.7734903665156936, 1.8658946539689207, 1.9562837130640909, 2.042751539057915, 2.122892653466968, 2.193686520478429], [0.0015120003493174115, 0.0016453382153351044, 0.0017900482821541844, 0.001947062962041602, 0.0021173877500562273, 0.0023021066886859875, 0.002502388226050207, 0.0027194914960287793, 0.002954773050894019, 0.00320969407943684, 0.0034858281462160916, 0.0037848694904418214, 0.004108641926153563, 0.0044591083888016044, 0.004838381177111012, 0.005248732943238343, 0.00569260848875382, 0.006172637428935729, 0.006691647793291092, 0.007252680636160278, 0.007859005737774902, 0.00851413848326635, 0.009221858014926577, 0.009986226761560834, 0.010811611458109232, 0.011702705778918946, 0.012664554719193612, 0.013702580871307051, 0.01482261275592704, 0.016030915382333394, 0.017334223228022094, 0.018739775844752785, 0.020255356316712284, 0.021889332816526143, 0.0236507035265444, 0.025549145216249236, 0.02759506579186571, 0.029799661161384567, 0.03217497678729623, 0.03473397433044651, 0.037490603821589835, 0.04045988183245383, 0.043657976155409714, 0.04710229754011029, 0.05081159907658642, 0.05480608385709701, 0.05910752159322928, 0.06373937490995818, 0.0687269360840703, 0.07409747503983642, 0.07988039945916516, 0.08610742790550369, 0.09281277689894973, 0.10003336291246825, 0.10780902028332434, 0.11618273604674655, 0.12520090269656434, 0.13491358985522606, 0.14537483578713253, 0.15664295960698899, 0.168780894909405, 0.181856545365408, 0.1959431625811127, 0.211119746175137, 0.22747146558150025, 0.24509010249516272, 0.2640745121123525, 0.2845311003329514, 0.306574312831917, 0.33032713030157135, 0.3559215621298035, 0.38349912820201587, 0.4132113152607084, 0.44521999015472774, 0.47969774714526947, 0.5168281599367066, 0.5568059009256507, 0.5998366798795802, 0.646136941320882, 0.6959332436106462, 0.7494612222194716, 0.806964013820291, 0.8686899852124572, 0.9348895698614175, 1.0058109626745562, 1.081694357526662, 1.1627643281288407, 1.2492198461090518, 1.3412212941770787, 1.4388736585944866, 1.5422048629454261, 1.6511379201610161, 1.7654552132406107, 1.8847527426512045, 2.008381567719428, 2.135372877836839, 2.2643421004643436, 2.3933661115491027, 2.51982585964305, 2.6402044131463347], [0.000311809967228322, 0.0003452371262278347, 0.0003821389868698809, 0.0004228633295298757, 0.0004677911428928001, 0.0005173396113769373, 0.0005719653513787946, 0.0006321679150657027, 0.000698493581677693, 0.0007715394576061729, 0.0008519579078977775, 0.0009404613432925704, 0.0010378273884515679, 0.001144904458665663, 0.0012626177740729782, 0.0013919758422523787, 0.0015340774420160254, 0.0016901191433033153, 0.0018614034002940817, 0.002049347257222868, 0.002255491708904247, 0.002481511760687407, 0.0027292272354668403, 0.003000614378506446, 0.003297818314211878, 0.0036231664126393478, 0.003979182627490375, 0.004368602871647783, 0.004794391501000359, 0.005259758982428298, 0.005768180827431334, 0.006323417879035601, 0.006929538046378621, 0.007590939588818913, 0.008312376059627773, 0.009098983028387077, 0.009956306711237019, 0.010890334649201581, 0.011907528587086201, 0.013014859719022371, 0.014219846481766643, 0.015530595093501163, 0.016955843054286268, 0.018505005844657956, 0.020188227081318136, 0.02201643241361972, 0.024001387471785297, 0.026155760207701443, 0.028493188001874534, 0.03102834994587326, 0.033777044748447026, 0.03675627475558232, 0.03998433662005812, 0.04348091920452015, 0.047267209353519024, 0.05136600622398556, 0.05580184491967411, 0.06060113023230214, 0.06579228134919164, 0.07140588844238602, 0.07747488210503245, 0.08403471664395898, 0.09112356826838254, 0.09878254922765005, 0.10705593893799598, 0.11599143308923839, 0.12564041162362163, 0.13605822631306738, 0.14730450840494919, 0.1594434964302602, 0.17254438373283484, 0.18668168453352876, 0.20193561632323254, 0.2183924949977675, 0.2361451372948256, 0.2552932626231887, 0.2759438830985461, 0.2982116662714667, 0.3222192493288489, 0.3480974760476651, 0.3759855179242467, 0.406030827964064, 0.43838885863402005, 0.4732224531808671, 0.5107007902216734, 0.5509977229788884, 0.5942893037836772, 0.6407502175254558, 0.6905487592452132, 0.7438398738762887, 0.8007556205387099, 0.8613922166651831, 0.925792540728495, 0.9939226020821446, 1.0656399890132333, 1.1406516355504273, 1.2184573402796646, 1.2982742380635885, 1.3789357448567068, 1.4587561938629154], [0.0001342725721263117, 0.00014603727525457771, 0.00015880576524531545, 0.00017266138310671186, 0.00018769421635013867, 0.00020400163477323104, 0.00022168886844931127, 0.00024086963126260122, 0.00026166679359513077, 0.0002842131080588479, 0.0003086519924767123, 0.00033513837465131124, 0.0003638396038205235, 0.00039493643408905173, 0.0004286240855442817, 0.0004651133892172801, 0.000504632022537223, 0.000547425842452816, 0.0005937603239602584, 0.0006439221123869805, 0.000698220698437356, 0.0007569902257142327, 0.0008205914411927007, 0.0008894137999442148, 0.0009638777362948396, 0.001044437114556288, 0.0011315818734983, 0.0012258408798421054, 0.0013277850072544277, 0.0014380304586173603, 0.0015572423507503339, 0.0016861385822758727, 0.0018254940069617168, 0.0019761449366502773, 0.0021389939998156496, 0.0023150153838836993, 0.002505260491728745, 0.002710864045240096, 0.0029330506715540246, 0.0031731420104954265, 0.0034325643849942127, 0.0037128570797643474, 0.004015681277389892, 0.004342829705189377, 0.004696237050866266, 0.005077991210043338, 0.005490345434369945, 0.0059357314550354855, 0.0064167736632759625, 0.006936304436882715, 0.007497380709876483, 0.00810330189146132, 0.008757629250186881, 0.00946420688998964, 0.010227184456513615, 0.011051041724879734, 0.011940615233918109, 0.012901127146810095, 0.013938216534087226, 0.015057973291928543, 0.016266974926545427, 0.017572326453901867, 0.01898170368271087, 0.02050340016701511, 0.02214637813188485, 0.023920323690702086, 0.025835706683531974, 0.02790384547099443, 0.030136977013795533, 0.03254833255051671, 0.03515221914977805, 0.037964107349903686, 0.04100072499950361, 0.044280157262188306, 0.047821952529457645, 0.05164723367277527, 0.055778813625521166, 0.06024131367287777, 0.06506128198211983, 0.07026730874539207, 0.07589013272125394, 0.08196273180019381, 0.08852038727946378, 0.09560070753693478, 0.10324359136389943, 0.11149110383470481, 0.12038722754451338, 0.1299774383529147, 0.1403080360698595, 0.1514251349194133, 0.16337318347344312, 0.17619283535345495, 0.18991792512593042, 0.2045712110543479, 0.22015841716027126, 0.23665992625381288, 0.25401922140400823, 0.27212681408152806, 0.29079788669773615, 0.3097411475137388]];

#[cfg(feature = "std")]
#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    static CUTOFF_TEST_INPUT: [f32; 101] = [0.0, 0.5877852522924732, -0.9510565162951536, 0.9510565162951535, -0.5877852522924728, -4.898587196589413e-16, 0.5877852522924736, -0.9510565162951538, 0.9510565162951533, -0.5877852522924725, -9.797174393178826e-16, 0.587785252292474, -0.951056516295154, 0.9510565162951532, -0.587785252292472, -1.4695761589768238e-15, 0.5877852522924744, -0.9510565162951541, 0.951056516295153, -0.5877852522924716, -1.959434878635765e-15, 0.5877852522924748, -0.9510565162951542, 0.9510565162951529, -0.5877852522924712, -2.4492935982947065e-15, 0.5877852522924751, -0.9510565162951544, 0.9510565162951528, -0.5877852522924708, -2.9391523179536475e-15, 0.5877852522924756, -0.9510565162951545, 0.9510565162951525, -0.5877852522924705, -3.429011037612589e-15, 0.587785252292476, -0.9510565162951546, 0.9510565162951524, -0.58778525229247, -3.91886975727153e-15, 0.5877852522924764, -0.9510565162951549, 0.9510565162951523, -0.5877852522924697, -4.408728476930472e-15, 0.5877852522924768, -0.951056516295155, 0.9510565162951521, -0.5877852522924693, -4.898587196589413e-15, 0.5877852522924887, -0.9510565162951552, 0.9510565162951563, -0.5877852522924688, -1.9599300631450357e-14, 0.5877852522924776, -0.9510565162951509, 0.9510565162951519, -0.5877852522924569, -5.878304635907295e-15, 0.5877852522924665, -0.9510565162951554, 0.9510565162951473, -0.587785252292468, 7.842691359635767e-15, 0.5877852522924784, -0.95105651629516, 0.9510565162951515, -0.5877852522924791, -6.858022075225178e-15, 0.5877852522924902, -0.9510565162951558, 0.9510565162951558, -0.5877852522924673, -2.1558735510086122e-14, 0.5877852522924791, -0.9510565162951515, 0.9510565162951512, -0.5877852522924554, -7.83773951454306e-15, 0.587785252292468, -0.9510565162951561, 0.9510565162951466, -0.5877852522924665, 5.883256481000002e-15, 0.5877852522924799, -0.9510565162951606, 0.9510565162951509, -0.5877852522924776, -8.817456953860943e-15, 0.5877852522924918, -0.9510565162951563, 0.9510565162951552, -0.5877852522924657, -2.3518170388721888e-14, 0.5877852522924807, -0.9510565162951521, 0.9510565162951506, -0.5877852522924538, -9.797174393178826e-15];
    static CUTOFF_TEST_OUTPUT: [f32; 101] = [0.0, 0.20817248253695878, 0.0846312532364365, -0.3097177262341411, 0.25898286035029705, 0.009372973898095042, -0.3071835272951865, 0.45106437907150543, -0.3557949885651552, 0.06422673819947913, 0.28448059836470707, -0.5254550446906006, 0.5437913795855674, -0.32378509991570786, -0.04635785092967576, 0.4134912277214284, -0.6242290513599532, 0.5882990116897372, -0.3152187207315345, -0.08964292716169187, 0.46724649233951626, -0.6680611990985434, 0.6111653921566373, -0.3161853514828968, -0.10416885760195049, 0.48782482797715243, -0.6861945217728828, 0.6217539122794384, -0.3181315022297383, -0.10883463947612498, 0.4955607848503048, -0.6935561523585955, 0.6264769606443668, -0.3195007641500944, -0.1102304559156358, 0.4984202394562113, -0.6965048623821712, 0.628534270362899, -0.32027567701321225, -0.11059418179821112, 0.4994549703401429, -0.6976698259794524, 0.6294128480519041, -0.32067190806800033, -0.11066009296140611, 0.4998196415871752, -0.6981235286771428, 0.6297815508463205, -0.3208625585477992, -0.11065438425729918, 0.49994386326912027, -0.6982975812762418, 0.6299338519808945, -0.3209505636599039, -0.11063935640788791, 0.49998423500234707, -0.6983632708114135, 0.6299958415761991, -0.3209899525352162, -0.11062835262340301, 0.4999964468500707, -0.6983876135858527, 0.630020717006998, -0.32100715583615674, -0.11062202380323702, 0.499999695163008, -0.6983964447020903, 0.6300305599250433, -0.3210145177544964, -0.11061875990307901, 0.5000003244858109, -0.698399566921172, 0.6300343993837549, -0.321017612838766, -0.11061718084157393, 0.5000003058095993, -0.6984006349567675, 0.6300358747977989, -0.3210188934600926, -0.11061644909099203, 0.5000001907297235, -0.6984009841637965, 0.6300364326737595, -0.321019415519255, -0.11061612059409717, 0.5000001024753468, -0.698401090822341, 0.6300366398454351, -0.32101962533455286, -0.11061597677110113, 0.5000000508476297, -0.6984011197350304, 0.6300367151918165, -0.32101970848359895, -0.11061591509625837, 0.5000000239816723, -0.6984011256611126, 0.6300367419130085, -0.32101974096917396, -0.11061588911943968, 0.5000000109097195];
    const STEP_TEST_MIN_OUTPUT: f32 = 1.0000000090010597;
    const STEP_TEST_MAX_OUTPUT: f32 = 1.0711971516552137;

    #[test]
    fn test() {
        let order = 5;
        println!("order {order}");
        let mut filter = butter5(0.4).unwrap();
        let out = (0..CUTOFF_TEST_INPUT.len()).map(|i| {filter.update(CUTOFF_TEST_INPUT[i])}).collect::<Vec<f32>>();
        // Check overall match to reference output to catch phase error, etc
        (0..CUTOFF_TEST_INPUT.len()).for_each(|i| { let expected = CUTOFF_TEST_OUTPUT[i]; let rel_err = (out[i] - expected).abs() / expected.abs().max(1e-4); assert!(rel_err < 0.05); });
        // Check approximate attenuation at cutoff frequency; should be -3dB or 1/sqrt(2) magnitude
        let maxmag = out.iter().fold(0.0_f32, |a, b| a.abs().max(b.abs()));
        let attenuation_rel_err = (maxmag - 0.707).abs() / 0.707;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);
        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter5(MIN_CUTOFF_RATIO).unwrap();
        (0..169).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        let step_min_rel_err = (step_min_final - STEP_TEST_MIN_OUTPUT).abs() / STEP_TEST_MIN_OUTPUT;
        println!("order {order} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);
        let mut filtermax = butter5(MAX_CUTOFF_RATIO).unwrap();
        (0..1).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - STEP_TEST_MAX_OUTPUT).abs() / STEP_TEST_MAX_OUTPUT;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }
}
