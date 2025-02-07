//! Butterworth filter of order 4.
//! Region of validity: cutoff ratio from 3.16e-02 to 4.00e-01 .
//! This file is autogenerated.
#![allow(clippy::style)]

use crate::SisoIirFilter;

/// Minimum tabulated cutoff ratio
#[allow(dead_code)]
pub const MIN_CUTOFF_RATIO: f64 = 0.03162277660168379;

/// Maximum tabulated cutoff ratio
#[allow(dead_code)]
pub const MAX_CUTOFF_RATIO: f64 = 0.4;

/// Initialise a Butterworth filter of order 4 by interpolating the coefficients from stored tables.
/// Cutoff ratio is the dimensionless ratio of the cutoff frequency to the sampling frequency.
/// Region of validity: cutoff ratio from 3.16e-02 to 4.00e-01
pub fn butter4(cutoff_ratio: f64) -> Result<SisoIirFilter<4>, &'static str> {
    let avals = &[&AVALS[0][..], &AVALS[1][..], &AVALS[2][..], &AVALS[3][..]];
    let cvals = &[&CVALS[0][..], &CVALS[1][..], &CVALS[2][..], &CVALS[3][..]];
    SisoIirFilter::new_interpolated(cutoff_ratio, &LOG10_CUTOFF_RATIOS, avals, cvals, &DVALS)
}

/// [dimensionless] Log base-10 of cutoff ratios, to improve float precision during interpolation
#[rustfmt::skip]
static LOG10_CUTOFF_RATIOS: [f64; 100] = [-1.5, -1.4888680808956771, -1.4777361617913543, -1.4666042426870314, -1.4554723235827085, -1.4443404044783856, -1.433208485374063, -1.42207656626974, -1.4109446471654172, -1.3998127280610944, -1.3886808089567715, -1.3775488898524486, -1.3664169707481257, -1.3552850516438029, -1.34415313253948, -1.3330212134351571, -1.3218892943308345, -1.3107573752265116, -1.2996254561221887, -1.2884935370178658, -1.277361617913543, -1.26622969880922, -1.2550977797048972, -1.2439658606005743, -1.2328339414962515, -1.2217020223919286, -1.2105701032876057, -1.199438184183283, -1.18830626507896, -1.1771743459746373, -1.1660424268703145, -1.1549105077659916, -1.1437785886616687, -1.1326466695573458, -1.121514750453023, -1.1103828313487, -1.0992509122443774, -1.0881189931400543, -1.0769870740357317, -1.0658551549314088, -1.054723235827086, -1.043591316722763, -1.0324593976184402, -1.0213274785141173, -1.0101955594097944, -0.9990636403054716, -0.9879317212011488, -0.9767998020968259, -0.9656678829925031, -0.9545359638881802, -0.9434040447838573, -0.9322721256795345, -0.9211402065752117, -0.9100082874708888, -0.8988763683665659, -0.887744449262243, -0.8766125301579202, -0.8654806110535974, -0.8543486919492744, -0.8432167728449517, -0.8320848537406288, -0.8209529346363059, -0.8098210155319832, -0.7986890964276603, -0.7875571773233374, -0.7764252582190145, -0.7652933391146917, -0.7541614200103689, -0.743029500906046, -0.7318975818017232, -0.7207656626974003, -0.7096337435930774, -0.6985018244887546, -0.6873699053844318, -0.6762379862801089, -0.665106067175786, -0.6539741480714631, -0.6428422289671403, -0.6317103098628175, -0.6205783907584946, -0.6094464716541718, -0.5983145525498489, -0.587182633445526, -0.5760507143412033, -0.5649187952368804, -0.5537868761325575, -0.5426549570282346, -0.5315230379239118, -0.5203911188195889, -0.5092591997152661, -0.49812728061094314, -0.4869953615066205, -0.4758634424022976, -0.46473152329797474, -0.45359960419365186, -0.442467685089329, -0.4313357659850061, -0.42020384688068324, -0.40907192777636037, -0.3979400086720376];

/// State-Space `D` 1x1 matrix
#[rustfmt::skip]
static DVALS: [f64; 100] = [7.609661935031952e-05, 8.380122457436599e-05, 9.227284867815485e-05, 0.00010158622362100626, 0.00011182316367792747, 0.00012307321290216068, 0.0001354343490409389, 0.0001490137485251852, 0.00016392861750121204, 0.000180307094243171, 0.00019828922868041815, 0.0002180280451956391, 0.0002396906952990662, 0.0002634597072631756, 0.00028953434031265873, 0.00031813205150825794, 0.0003494900840423613, 0.00038386718628145204, 0.0004215454715481422, 0.00046283242933650937, 0.0005080630994018316, 0.0005576024209632275, 0.0006118477701089365, 0.000671231699403506, 0.0007362248946688966, 0.0008073393649531053, 0.0008851318828166409, 0.0009702076932663859, 0.0010632245109562048, 0.0011648968266635728, 0.0012760005455524095, 0.001397377981356278, 0.0015299432323778498, 0.0016746879671158742, 0.0018326876494189493, 0.002005108235347669, 0.0021932133764279815, 0.002398372166727619, 0.0026220674742169128, 0.0028659049002235917, 0.003131622414500844, 0.003421100717550759, 0.003736374386437162, 0.004079643865450943, 0.0044532883687324306, 0.004859879768397982, 0.005302197548962591, 0.0057832449170137955, 0.0063062661643086305, 0.006874765392889501, 0.007492526722625371, 0.008163636114986628, 0.008892504962093787, 0.009683895607415503, 0.010542948984248433, 0.0114752145806574, 0.012486682965316293, 0.013583821138162554, 0.014773611003538483, 0.01606359130221013, 0.01746190338311603, 0.01897734124681542, 0.020619406351458345, 0.02239836773994345, 0.024325328125252548, 0.026412296661501702, 0.028672269233097843, 0.031119317215991402, 0.033768685806282756, 0.03663690317584738, 0.039741901906335565, 0.04310315437684484, 0.04674182404271471, 0.050680934850435164, 0.05494556139525205, 0.059563042854190035, 0.06456322423063425, 0.0699787290428383, 0.07584526829671197, 0.08220199142623888, 0.08909188589144809, 0.0965622333292968, 0.1046651316007602, 0.11345809382215596, 0.12300473757793845, 0.1333755800706802, 0.1446489580783449, 0.15691209539396309, 0.17026234508938567, 0.18480863969095487, 0.20067318946020296, 0.21799347779745373, 0.2369246137979589, 0.25764211579925744, 0.2803452171674747, 0.3052608076427798, 0.33264815172150725, 0.3628045617112942, 0.3960722498666854, 0.43284664499029174];

/// State-Space `A` matrix, first row
#[rustfmt::skip]
static AVALS: [[f64; 100]; 4] = [[3.48120132755602, 3.4677531427727724, 3.453957488009553, 3.4398054718794926, 3.4252879818821533, 3.4103956793557586, 3.3951189943430995, 3.3794481203711584, 3.363373009144521, 3.346883365152678, 3.3299686401913187, 3.312618027797698, 3.294820457600147, 3.276564589581704, 3.2578388082577563, 3.2386312167674864, 3.218929630878666, 3.1987215729051903, 3.177994265536425, 3.1567346255770747, 3.134929257595888, 3.1125644474809775, 3.0896261558989444, 3.066100011654245, 3.041971304944446, 3.0172249805059694, 2.991845630643876, 2.965817488137851, 2.93912441901512, 2.9117499151793025, 2.8836770868822614, 2.854888655023844, 2.825366943261926, 2.795093869912391, 2.7640509396156023, 2.7322192347424092, 2.6995794065089154, 2.6661116657649013, 2.6317957734161066, 2.596611030435258, 2.5605362674110443, 2.523549833577838, 2.485629585262054, 2.446752873673449, 2.40689653196141, 2.3660368614472906, 2.324149616934161, 2.2812099909847734, 2.237192597047292, 2.1920714512960706, 2.145819953041804, 2.098410863551382, 2.0498162831029694, 2.0000076260860933, 1.9489555939398628, 1.8966301457049544, 1.8430004659466235, 1.7880349297868814, 1.731701064764126, 1.6739655092181152, 1.6147939668772768, 1.5541511573043092, 1.492000761834892, 1.4283053646236277, 1.3630263883913558, 1.2961240244492434, 1.2275571565583383, 1.157283278169216, 1.0852584025761334, 1.0114369655149582, 0.935771719735599, 0.8582136210898827, 0.7787117056971056, 0.6972129577850369, 0.6136621678577932, 0.5280017809184376, 0.4401717345793863, 0.3501092870349292, 0.2577488350565878, 0.16302172241460366, 0.06585603944153627, -0.03382358514661334, -0.13609619854931232, -0.24104474203849716, -0.34875626990925435, -0.4593221651579887, -0.5728383465775577, -0.6894054604212425, -0.8091290478983957, -0.9321196774590252, -1.05849302803266, -1.1883699060385355, -1.3218761750185477, -1.4591425721235503, -1.600304380420097, -1.7455009201805858, -1.8948748162334104, -2.048570992595213, -2.2067353409053134, -2.3695130071820376], [-4.5739752781920595, -4.540364185021797, -4.506054818140668, -4.471037794556154, -4.435303912639332, -4.398844176328249, -4.3616498207581635, -4.323712339380616, -4.285023512635141, -4.245575438239211, -4.205360563163781, -4.164371717363656, -4.122602149333634, -4.0800455635632416, -4.036696159964599, -3.992548675349927, -3.9475984270368127, -3.901841358661417, -3.8552740882815733, -3.807893958853661, -3.759699091169184, -3.710688439339034, -3.6608618489156, -3.610220117745235, -3.558765059646092, -3.506499571009081, -3.4534277004227385, -3.3995547214261137, -3.3448872084975982, -3.289433116391899, -3.233201862942195, -3.176204415450203, -3.1184533807932597, -3.0599630993850147, -3.0007497431350303, -2.9408314175625123, -2.880228268231143, -2.8189625916853975, -2.7570589510845247, -2.6945442967483864, -2.6314480918505554, -2.5678024435182847, -2.5036422396271143, -2.4390052916103073, -2.3739324836405693, -2.3084679285844447, -2.2426591311789834, -2.1765571589366624, -2.110216821349149, -2.043696858034113, -1.977060136553524, -1.9103738607276215, -1.8437097903778263, -1.7771444735558009, -1.7107594924565708, -1.644641724373236, -1.578883619231665, -1.5135834954483969, -1.4488458560867477, -1.3847817275482732, -1.3215090233330584, -1.2591529357372453, -1.197846358734539, -1.1377303457158285, -1.0789546062438713, -1.0216780465251944, -0.966069358917433, -0.9123076664861662, -0.8605832294117113, -0.8110982209350798, -0.7640675815369793, -0.7197199611797381, -0.6782987607266149, -0.6400632851057336, -0.6052900224287925, -0.5742740651322102, -0.5473306913076768, -0.5247971267599671, -0.5070345110046985, -0.4944300934318514, -0.48739968924807553, -0.48639042860784437, -0.49188383658395385, -0.5043992863400362, -0.5244978730696702, -0.5527867619590641, -0.5899240695861752, -0.6366243447191505, -0.6936647212884226, -0.7618918231544796, -0.8422295068187515, -0.9356875338804884, -1.0433712690149153, -1.166492500355618, -1.3063814757219898, -1.4645002377745677, -1.642457320640145, -1.8420238353914853, -2.0651509161026143, -2.31398841441588], [2.68595145033013, 2.657659118658445, 2.6289053742418926, 2.599690148435927, 2.570013821594103, 2.5398772464716637, 2.5092817721027547, 2.4782292681156677, 2.4467221494454963, 2.414763401398238, 2.3823566050146745, 2.349505962676319, 2.3162163238892584, 2.2824932111749745, 2.248342845989924, 2.2137721745882537, 2.1787888937338202, 2.143401476159472, 2.1076191956626698, 2.071452151717259, 2.034911293471636, 1.998008442993436, 1.9607563176103595, 1.9231685511858303, 1.8852597141567942, 1.8470453321491358, 1.808541902973937, 1.7697669117950399, 1.730738844245166, 1.6914771972541753, 1.6520024873387587, 1.6123362560881769, 1.572501072565233, 1.5325205323256945, 1.492419252742725, 1.4522228643053166, 1.4119579975414096, 1.3716522651969167, 1.3313342392813354, 1.2910334225685716, 1.2507802141180893, 1.2106058683558978, 1.1705424472271424, 1.130622764901583, 1.0908803244795522, 1.0513492461085487, 1.0120641858786161, 0.9730602448173367, 0.9343728672515985, 0.8960377277420312, 0.8580906057258916, 0.8205672469234327, 0.7835032104695916, 0.7469337006248543, 0.7108933817937828, 0.675416175433773, 0.6405350372664912, 0.6062817130058311, 0.5726864705840186, 0.5397778065857759, 0.5075821242821432, 0.4761233802824414, 0.4454226963850604, 0.41549793269376856, 0.3863632174623014, 0.358028428419746, 0.3304986194933913, 0.30377338586098845, 0.27784615910319477, 0.2527034228561651, 0.22832383774401208, 0.20467726245332032, 0.1817236555392251, 0.15941183985490937, 0.13767810828898752, 0.11644464567514493, 0.09561773718016539, 0.07508572802689387, 0.05471669288025611, 0.034355765387101646, 0.013822068931467306, -0.0070948217028244, -0.028638971811124122, -0.051092832082258816, -0.07478275865173366, -0.10008545339279376, -0.1274354999463237, -0.15733420767078557, -0.1903600206370373, -0.22718080397696677, -0.2685683877848747, -0.315415832438249, -0.3687579824703646, -0.4297960037265253, -0.4999267563096064, -0.5807780509161227, -0.6742510772259451, -0.7825715903061801, -0.9083518063819276, -1.0546654058785674], [-0.5943950456036955, -0.5863888960026092, -0.5782844096896281, -0.5700832053372001, -0.5617870614557704, -0.5533979209056072, -0.5449178952723458, -0.5363492690826133, -0.5276945038348946, -0.5189562418195957, -0.510137309701099, -0.501240721833491, -0.49226968328055637, -0.48322759250964753, -0.47411804372808386, -0.46494482882994587, -0.455711938920351, -0.4464235653837483, -0.43708410046229157, -0.4276981373100571, -0.41827046948876934, -0.4088060898707917, -0.3993101889154465, -0.3897881522852968, -0.3802455577698504, -0.3706881714852733, -0.3611219433201406, -0.3515530015990391, -0.34198764693798706, -0.3324323452681963, -0.3228937200076639, -0.31337854336351856, -0.30389372675194487, -0.2944463103269249, -0.28504345161400024, -0.27569241325077576, -0.26640054984202965, -0.2571752939440628, -0.24802414120038835, -0.23895463465902073, -0.2299743483105912, -0.22109086989626334, -0.2123117830450767, -0.20364464881194003, -0.19509698670011177, -0.18667625526576237, -0.1783898324171952, -0.1702449955376683, -0.1622489015786796, -0.1544085672902212, -0.14673084977617765, -0.13922242758697895, -0.131889782588235, -0.12473918287379454, -0.1177766670250497, -0.11100803005600957, -0.10443881142651029, -0.09807428555491654, -0.09191945531801271, -0.08597904909097978, -0.08025752195621819, -0.0747590617985518, -0.06948760110874685, -0.06444683544066293, -0.05964024961382673, -0.05507115292782217, -0.05074272486386257, -0.046658072999900484, -0.042820305168141164, -0.03923261824960142, -0.03589840644400089, -0.03282139239298246, -0.03000578519315076, -0.027456470141174916, -0.025179236042020985, -0.023181047128412798, -0.02147036814202249, -0.020057552987268433, -0.018955309679537145, -0.018179257189676497, -0.017748593388097825, -0.017686897811466873, -0.018023098667773, -0.018792640693702463, -0.020038899616357467, -0.02181490062103652, -0.02418541314346163, -0.02722951349223044, -0.031043731606315544, -0.03574593046480558, -0.041480108726961076, -0.04842237240198769, -0.056788394263515785, -0.0668427765824253, -0.07891086222790158, -0.09339371341320105, -0.11078721344461678, -0.13170656908783035, -0.15691793447711186, -0.18737949236818488]];

/// State-Space `C` vector
#[rustfmt::skip]
static CVALS: [[f64; 100]; 4] = [[0.0005692941297061355, 0.0006258068581834284, 0.0006877978913445047, 0.0007557817423631361, 0.0008303191933537558, 0.0009120212051346029, 0.0010015531270791403, 0.0010996392264636383, 0.001207067557534801, 0.0013246951913141692, 0.00145345382791519, 0.0015943558138631219, 0.0017485005875640318, 0.0019170815766527852, 0.002101393571444523, 0.0023028405991019577, 0.002522944323391677, 0.0027633529950147045, 0.003025850977435413, 0.0033123688728725516, 0.003624994272626976, 0.003965983155172374, 0.00433777195439276, 0.004742990318977813, 0.005184474582244117, 0.005665281959494756, 0.006188705487415118, 0.006758289716880875, 0.00737784716687161, 0.008051475542884587, 0.008783575718268386, 0.009578870471179269, 0.010442423963339074, 0.01137966193936512, 0.012396392617074155, 0.013498828229747987, 0.01469360717079679, 0.015987816679468823, 0.01738901599312357, 0.018905259876993402, 0.02054512242617013, 0.02231772101663128, 0.02423274026229221, 0.026300455813160123, 0.028531757805455894, 0.030938173747823473, 0.03353189059818102, 0.03632577575305896, 0.039333396635035614, 0.04256903852366929, 0.04604772023060992, 0.04978520716971485, 0.05379802131724895, 0.05810344749471464, 0.0627195353364673, 0.06766509622473746, 0.07295969438447086, 0.07862363122766225, 0.08467792191939252, 0.09114426300291639, 0.0980449897657138, 0.10540302184855867, 0.11324179539079245, 0.12158517976154783, 0.13045737664200785, 0.13988279888985972, 0.1498859262242481, 0.16049113430607587, 0.1717224932403528, 0.18360353087743392, 0.19615695551779738, 0.20940433170552758, 0.22336570169855519, 0.23805914389212332, 0.2535002579009821, 0.2697015641206937, 0.2866718033221731, 0.30441511910415187, 0.3229301027348797, 0.34220867593317145, 0.3622347823169804, 0.3829828523962266, 0.40441599987151383, 0.4264838983310827, 0.4491202768529053, 0.4722399600654533, 0.4957353623336141, 0.5194723262051153, 0.54328517118242, 0.5669707891434455, 0.59028158588411, 0.6129170224626412, 0.634513452936843, 0.6546318836623477, 0.6727431896069656, 0.6882102099355782, 0.7002660015223543, 0.7079873457421997, 0.7102623681340471, 0.7057508245415649], [0.00010851566043956415, 0.0001223192687277761, 0.00013785057769916473, 0.0001553214965102855, 0.00017496926668314461, 0.00019705939157628786, 0.0002218888900267237, 0.00024978990791540783, 0.000281133724620806, 0.00031633519479993665, 0.00035585766968971374, 0.00040021844616903025, 0.00044999479617919376, 0.0005058306337822637, 0.0005684438821579596, 0.0006386346082138976, 0.0007172939982235788, 0.0008054142540227561, 0.0009040994957968634, 0.0010145777643869984, 0.0011382142233333245, 0.0012765256685636594, 0.0014311964617177305, 0.0016040960115661766, 0.001797297936824086, 0.002013101052851824, 0.0022540523342535245, 0.0025229720151906344, 0.0028229809992787046, 0.0031575307611744464, 0.00353043593231915, 0.003945909773700967, 0.004408602748836624, 0.004923644420336565, 0.005496688903273041, 0.006133964117962262, 0.006842325093517346, 0.007629311581421102, 0.008503210245164147, 0.009473121697420817, 0.010549032659970389, 0.01174189352325556, 0.01306370158167805, 0.014527590216985279, 0.016147924294841947, 0.017940402028264753, 0.019922163545279765, 0.022111906376072334, 0.024530008046022733, 0.027198655924167082, 0.030141984430387583, 0.03338621964735733, 0.036959831312967, 0.04089369208328326, 0.045221243852202234, 0.049978670788459106, 0.055205078599420795, 0.06094267934912942, 0.0672369809393166, 0.07413698009915681, 0.08169535741333829, 0.08996867253747737, 0.09901755729138788, 0.1089069037674246, 0.11970604392238045, 0.13148891631164325, 0.1443342146418602, 0.15832551162398453, 0.1735513501535363, 0.19010529206858368, 0.2080859125627597, 0.227596725666242, 0.24874602293401338, 0.27164660345001174, 0.2964153682822176, 0.3231727463736205, 0.35204191123259987, 0.38314773832103405, 0.41661544125742833, 0.4525688102562737, 0.4911279578506717, 0.5324064539194134, 0.5765077031162149, 0.6235203813795394, 0.6735127022305091, 0.726525225392398, 0.7825618464590922, 0.8415785124550879, 0.9034690883839722, 0.9680476467168884, 1.035026255370405, 1.103987087142392, 1.1743473478285082, 1.2453150989399573, 1.3158335044900007, 1.384510320480571, 1.4495285183366404, 1.5085327200068002, 1.5584845295451035, 1.595477748215432], [0.0005087783025104887, 0.0005579199869422719, 0.0006116679824992296, 0.0006704375992483641, 0.0007346797309383631, 0.0008048837047090034, 0.0008815803395287858, 0.0009653452270474826, 0.0010568022493730423, 0.0011566273491635563, 0.001265552568371732, 0.0013843703730003726, 0.0015139382823323284, 0.0016551838222990467, 0.0018091098239610133, 0.00197680008950669, 0.0021594254497510415, 0.0023582502388506555, 0.002574639213872105, 0.0028100649489796762, 0.003066115736376316, 0.0033445040287710184, 0.003647075461092654, 0.0039758184924658675, 0.004332874713154186, 0.004720549865309302, 0.005141325630998676, 0.005597872246177384, 0.006093062005090293, 0.006629983726109474, 0.007211958257307831, 0.007842555108225147, 0.008525610303389487, 0.00926524556330735, 0.010065888929932447, 0.010932296966169448, 0.011869578672874213, 0.012883221282187652, 0.013979118102998556, 0.015163598612986072, 0.01644346101214975, 0.017826007475106556, 0.019319082363805615, 0.0209311136887737, 0.022671158135613573, 0.024548950004275333, 0.026574954441608782, 0.028760425382843242, 0.031117468654831316, 0.033659110732962266, 0.03639937368433652, 0.0393533568717058, 0.04253732603529241, 0.04596881041217364, 0.04966670859448542, 0.053651403866979086, 0.0579448897997889, 0.06257090690146022, 0.06755509115755161, 0.07292513528783785, 0.0787109635456757, 0.08494492085046884, 0.0916619769807592, 0.0988999464514351, 0.10669972454130898, 0.11510553971068037, 0.12416522233167304, 0.13393048922034947, 0.14445724287436926, 0.15580588353877606, 0.16804163118784288, 0.18123485314833482, 0.1954613913024722, 0.2108028804718153, 0.22734704653278293, 0.24518796883724978, 0.2644262863285259, 0.2851693199879315, 0.3075310754386597, 0.33163207803674793, 0.3575989777538184, 0.38556384148848927, 0.41566302464951904, 0.4480354799525953, 0.482820316708443, 0.5201533648798199, 0.5601624200239509, 0.6029607413730804, 0.6486382368326169, 0.6972495834169388, 0.7487982828758364, 0.8032153169242221, 0.86033061261013, 0.9198349114348623, 0.9812287936044509, 1.043754453687264, 1.1063042322505832, 1.1672977040164327, 1.2245160558425836, 1.2748781974393046], [3.0865165821299295e-05, 3.466111701253679e-05, 3.891289884992769e-05, 4.367362364104142e-05, 4.900235715261696e-05, 5.4964752762931806e-05, 6.163374861397022e-05, 6.909033342044177e-05, 7.742438702456974e-05, 8.673560224132333e-05, 9.713449501868362e-05, 0.00010874351044183196, 0.00012169823263889854, 0.00013614870719909475, 0.00015226088529151968, 0.00017021819927443143, 0.0001902232802099804, 0.00021249982834785876, 0.00023729464831256995, 0.00026487986142259567, 0.00029555530828510824, 0.00032965115554676333, 0.00036753072143924236, 0.00040959353553769365, 0.00045627864895147265, 0.000508068211990557, 0.0005654913371992806, 0.0006291282665241081, 0.0006996148622875004, 0.0007776474425803217, 0.0008639879826671834, 0.0009594697050305935, 0.0010650030817716282, 0.0011815822742497063, 0.0013102920360982232, 0.0014523151071156658, 0.001608940127026673, 0.0017815700997621847, 0.0019717314407547916, 0.002181083641823166, 0.0024114295905711724, 0.002664726583904731, 0.0029430980783287336, 0.0032488462231934, 0.003584465227086077, 0.003952655612191606, 0.004356339416760292, 0.004798676411923539, 0.005283081406086762, 0.005813242718117039, 0.006393141909643831, 0.0070270748771214565, 0.007719674415978437, 0.008475934382311364, 0.009301235592268519, 0.010201373615588625, 0.011182588637759006, 0.012251597584931488, 0.013415628727013027, 0.01468245899705997, 0.016060454288948235, 0.017558613019772542, 0.019186613267809047, 0.02095486382006786, 0.022874559483924247, 0.024957741032881134, 0.02721736016418017, 0.029667349841620616, 0.03232270037493065, 0.03519954153970175, 0.038315230958844315, 0.04168844883366712, 0.04533929891095296, 0.04928941527598736, 0.05356207413541964, 0.05818230915067538, 0.06317702803796658, 0.06857512697707988, 0.0744075977484202, 0.0807076202825977, 0.08751063023458197, 0.0948543469759545, 0.10277874160684425, 0.1113259166311638, 0.12053985798927774, 0.1304660050461652, 0.14115056326644887, 0.15263945537533902, 0.1649767665457689, 0.1782024829072667, 0.19234924374280768, 0.207437716434341, 0.2234700454188692, 0.24042060141466434, 0.2582229343593209, 0.2767513672575077, 0.2957949899347794, 0.3150208176388856, 0.33392141051390256, 0.35174006037873895]];

#[cfg(feature = "std")]
#[cfg(test)]
#[rustfmt::skip]
mod test {
    use super::*;
    static CUTOFF_TEST_INPUT: [f32; 101] = [0.0, 0.5877852522924732, -0.9510565162951536, 0.9510565162951535, -0.5877852522924728, -4.898587196589413e-16, 0.5877852522924736, -0.9510565162951538, 0.9510565162951533, -0.5877852522924725, -9.797174393178826e-16, 0.587785252292474, -0.951056516295154, 0.9510565162951532, -0.587785252292472, -1.4695761589768238e-15, 0.5877852522924744, -0.9510565162951541, 0.951056516295153, -0.5877852522924716, -1.959434878635765e-15, 0.5877852522924748, -0.9510565162951542, 0.9510565162951529, -0.5877852522924712, -2.4492935982947065e-15, 0.5877852522924751, -0.9510565162951544, 0.9510565162951528, -0.5877852522924708, -2.9391523179536475e-15, 0.5877852522924756, -0.9510565162951545, 0.9510565162951525, -0.5877852522924705, -3.429011037612589e-15, 0.587785252292476, -0.9510565162951546, 0.9510565162951524, -0.58778525229247, -3.91886975727153e-15, 0.5877852522924764, -0.9510565162951549, 0.9510565162951523, -0.5877852522924697, -4.408728476930472e-15, 0.5877852522924768, -0.951056516295155, 0.9510565162951521, -0.5877852522924693, -4.898587196589413e-15, 0.5877852522924887, -0.9510565162951552, 0.9510565162951563, -0.5877852522924688, -1.9599300631450357e-14, 0.5877852522924776, -0.9510565162951509, 0.9510565162951519, -0.5877852522924569, -5.878304635907295e-15, 0.5877852522924665, -0.9510565162951554, 0.9510565162951473, -0.587785252292468, 7.842691359635767e-15, 0.5877852522924784, -0.95105651629516, 0.9510565162951515, -0.5877852522924791, -6.858022075225178e-15, 0.5877852522924902, -0.9510565162951558, 0.9510565162951558, -0.5877852522924673, -2.1558735510086122e-14, 0.5877852522924791, -0.9510565162951515, 0.9510565162951512, -0.5877852522924554, -7.83773951454306e-15, 0.587785252292468, -0.9510565162951561, 0.9510565162951466, -0.5877852522924665, 5.883256481000002e-15, 0.5877852522924799, -0.9510565162951606, 0.9510565162951509, -0.5877852522924776, -8.817456953860943e-15, 0.5877852522924918, -0.9510565162951563, 0.9510565162951552, -0.5877852522924657, -2.3518170388721888e-14, 0.5877852522924807, -0.9510565162951521, 0.9510565162951506, -0.5877852522924538, -9.797174393178826e-15];
    static CUTOFF_TEST_OUTPUT: [f32; 101] = [0.0, 0.25442087442956923, 0.0031683041842727744, -0.304693914037045, 0.38625525728967247, -0.20113280307219114, -0.1335664473624022, 0.4316297436541429, -0.5368100794474008, 0.3916876922012213, -0.05607350125670594, -0.3253369384980502, 0.5881312810635955, -0.6178947675900981, 0.39674873556653684, -0.009559580584369535, -0.39107143736244776, 0.6461132579465025, -0.6531799059581529, 0.40679985421929776, -0.0005607237386393171, -0.4093031364389759, 0.6645013997494591, -0.6659370632184979, 0.41202029826642633, 0.0006121897697404734, -0.414153445533343, 0.6701647698795844, -0.6703481688124726, 0.4142631707516016, 0.00044080938088747923, -0.41534629736692524, 0.6718488907562263, -0.6718179891526779, 0.41513906570239184, 0.00021384408371079368, -0.4156007119607376, 0.6723285625990629, -0.672290617704896, 0.41545986103040733, 8.915115411642042e-05, -0.4156389827971401, 0.6724580130015309, -0.6724373949917936, 0.4155718328143772, 3.411780557658467e-05, -0.4156371373966426, 0.6724903971184946, -0.6724813392606996, 0.4156093887005343, 1.2314975904026604e-05, -0.4156321008127822, 0.6724975397236945, -0.6724939608659739, 0.41562154134513896, 4.248563680374351e-06, -0.4156291352302292, 0.6724987263777754, -0.6724974047724283, 0.4156253396053624, 1.4104839639685943e-06, -0.4156277896791777, 0.6724987447844556, -0.6724982804684272, 0.41562648467732777, 4.5202983917931127e-07, -0.4156272481756989, 0.6724986359563698, -0.6724984793461899, 0.41562681623738545, 1.399053674314013e-07, -0.415627045677226, 0.6724985659666731, -0.6724985150244858, 0.41562690765053745, 4.1731563088307747e-08, -0.4156269738387863, 0.6724985331926592, -0.6724985171878297, 0.41562693124552097, 1.1933478365780177e-08, -0.4156269494075323, 0.6724985197754942, -0.6724985149235411, 0.4156269367441324, 3.2375771993624305e-09, -0.41562694140016926, 0.6724985147007234, -0.6724985132867624, 0.4156269377930943, 8.163191925828859e-10, -0.4156269388656444, 0.6724985128848595, -0.6724985124918292, 0.41562693789204563, 1.8263838317966348e-10, -0.41562693809123, 0.6724985122628511, -0.6724985121602258, 0.41562693784759613, 3.1550962221039575e-11];
    const STEP_TEST_MIN_OUTPUT: f32 = 1.0000000000269045;
    const STEP_TEST_MAX_OUTPUT: f32 = 1.1385974695318566;

    #[test]
    fn test() {
        let order = 4;
        println!("order {order}");
        let mut filter = butter4(0.4).unwrap();
        let out = (0..CUTOFF_TEST_INPUT.len()).map(|i| {filter.update(CUTOFF_TEST_INPUT[i])}).collect::<Vec<f32>>();
        // Check overall match to reference output to catch phase error, etc
        (0..CUTOFF_TEST_INPUT.len()).for_each(|i| { let expected = CUTOFF_TEST_OUTPUT[i]; let rel_err = (out[i] - expected).abs() / expected.abs().max(1e-4); assert!(rel_err < 0.05); });
        // Check approximate attenuation at cutoff frequency; should be -3dB or 1/sqrt(2) magnitude
        let maxmag = out.iter().fold(0.0_f32, |a, b| a.abs().max(b.abs()));
        let attenuation_rel_err = (maxmag - 0.707).abs() / 0.707;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);
        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter4(MIN_CUTOFF_RATIO).unwrap();
        (0..309).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        let step_min_rel_err = (step_min_final - STEP_TEST_MIN_OUTPUT).abs() / STEP_TEST_MIN_OUTPUT;
        println!("order {order} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);
        let mut filtermax = butter4(MAX_CUTOFF_RATIO).unwrap();
        (0..1).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - STEP_TEST_MAX_OUTPUT).abs() / STEP_TEST_MAX_OUTPUT;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }
}
