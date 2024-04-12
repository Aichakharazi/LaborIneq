function [de, Sfe, ej, alp, N_forcef, e_m, u_m, ms, thigtness, Vs, zz, eeta, bbeta, aalpha, ssp, ssigmaparam, xiparam, ssigmaTFP, rrhoTFP, theta_hf, theta_hn, theta_lf, theta_ln, phi_hf, phi_hn, phi_lf, phi_ln, varrho, upsilon, varepsilon, zetaparam, ddeltaCov, zzV, rsV, oy_hf_initV, oy_hn_initV, oy_lf_initV, oy_ln_initV, om_hf_initV, om_hn_initV, om_lf_initV, om_ln_initV, LYV, LMV, l_hftmV, l_hntmV, l_lftmV, l_lntmV, l_hftyV, l_hntyV, l_lftyV, l_lntyV, yV, margcostV, lambdacV, cmV, cyV, coldV, e2V, omegamV, omegayV, SV, N_forceV, e_mV, u_mV, mV, thigtnessV, VsV, zzVV, rsVV, oy_hf_initVV, oy_hn_initVV, oy_lf_initVV, oy_ln_initVV, om_hf_initVV, om_hn_initVV, om_lf_initVV, om_ln_initVV, LYVV, LMVV, l_hftmVV, l_hntmVV, l_lftmVV, l_lntmVV, l_hftyVV, l_hntyVV, l_lftyVV, l_lntyVV, yVV, margcostVV, lambdacVV, cmVV, cyVV, coldVV, e2VV, omegamVV, omegayVV, SVV, N_forceVV, e_mVV, u_mVV, mVV, thigtnessVV, VsVV] = sstatelabor;
 
%%


e2V = 1;
e2VV = 1;


eeta       = 1.5;          % 1.5the curvature on the disutility of labor  1-4
bbeta      = 0.95;       % disocunt factor 
aalpha      = 0.42;       % CPS youth labor share
ssigmaTFP  = 0.08 ;        % standard deviation of TFP shock
rrhoTFP    = 0.98;       % quarterly autocorrelation of TFP shock
zz         = 1;
zzV         = 1;
zzVV         = 1;

inzz = 1;
ddeltaCov = -0.0017;
ssp_elas          =  0.076580  ;               % CPS estimated  (robust and with County and year FE)
ssp               = (1/(1+ssp_elas)) ;       % CPS
sigmaparam_elas   = 0.03924;                 % CPS estimated  (!!!! not robust and with County and year FE)
ssigmaparam       = (1/(sigmaparam_elas +1)); % CPS
xiparam_elas  = 0.07146;                       % CPS estimated  (robust and with County and year FE)
xiparam =(1/(xiparam_elas+1));               % CPS
varrho_elas      = 0.14092;                     % CPS estimated (robust  with and without County and year FE)
varrho           = 1/(varrho_elas +1);          % CPS 
upsilon_elas     = 0.11723 ;                    % CPS estimated (!!!!!not robust and with County and year FE)
upsilon          = 1/(upsilon_elas +1);         % CPS 
varepsilon_elas  = 0.05898;                      % CPS estimated (!!!!!not robust and with County and year FE)
varepsilon       = 1/(varepsilon_elas +1)  ;    % CPS 
zetaparam_elas   = 0.11216;                      % CPS estimated (not robust (significant only) without County and year FE)
zetaparam        = 1/(zetaparam_elas +1)  ;     % CPS
rs               = (1/bbeta)-1;
rsV               = (1/bbeta)-1;
rsVV               = (1/bbeta)-1;

zma = 0.5;
zya = 0.5;

e_m = 0.95;   % employmnent rate
e_mV = 0.95; 
e_mVV = 0.95; 

u_m = 1 - e_m;
u_mV = 1 - e_m;
u_mVV = 1 - e_m;


alp = 0.5;
de = 0.025;


theta_hf  = 0.351;    % CPS high skill foreign labor share
theta_hn  = 0.256;    % CPS high skill native labor share
theta_lf  = 0.242;    % CPS low skill foreign labor share
theta_ln  = 0.249;    % CPS low skill native labor share


phi_hf    = 0.263;      % CPS high skill foreign labor share
phi_hn    = 0.265;      % CPS high skill native labor share
phi_lf    = 0.238;      % CPS low skill foreign labor share
phi_ln    = 0.232;      % CPS low skill native labor share

N_forcef =1.428100000000000;

%lambdacf = 0.040120000000000;





%% 
%% seatdy states
lambdac =    0.0734;
cm =   12.9485;
cy =   12.9485;
cold =   12.9485;
LM =    0.7156;
LY =    0.6294;
y =    0.6795;
omegam =    8.6453;
omegay =    6.2599;
oy_hf_init =    2.1589;
oy_hn_init =    1.2628;
oy_lf_init =    1.5749;
oy_ln_init =    1.4290;
om_hf_init =    2.9619;
om_hn_init =    1.3332;
om_lf_init =    2.8993;
om_ln_init =    1.4348;
N_force =    1.4243;
margcost =   14.8932;
l_hfty =    0.2014;
l_hnty =    0.1259;
l_lfty =    0.1511;
l_lnty =    0.1385;
l_hftm =    0.2433;
l_hntm =    0.1145;
l_lftm =    0.2290;
l_lntm =    0.1217;
Sfe =    0.0388;
ms =    0.0553;
thigtness =    2.0288;
Vs =    0.0787;
ej =    0.9178;
%y_pi =   28.5354;


lambdacV =    0.0734;
cmV =   12.9485;
cyV =   12.9485;
coldV =   12.9485;
LMV =    0.7156;
LYV =    0.6294;
yV =    0.6795;
omegamV =    8.6453;
omegayV =    6.2599;
oy_hf_initV =    2.1589;
oy_hn_initV =    1.2628;
oy_lf_initV =    1.5749;
oy_ln_initV =    1.4290;
om_hf_initV =    2.9619;
om_hn_initV =    1.3332;
om_lf_initV =    2.8993;
om_ln_initV =    1.4348;
N_forceV =    1.4243;
margcostV =   14.8932;
l_hftyV =    0.2014;
l_hntyV =    0.1259;
l_lftyV =    0.1511;
l_lntyV =    0.1385;
l_hftmV =    0.2433;
l_hntmV =    0.1145;
l_lftmV =    0.2290;
l_lntmV =    0.1217;
SV =    0.0388;
mV =    0.0553;
thigtnessV =    2.0288;
VsV =    0.0787;
%ejV =    0.9178;
%y_piV =   28.5354;




lambdacVV =    0.0734;
cmVV =   12.9485;
cyVV =   12.9485;
coldVV =   12.9485;
LMVV =    0.7156;
LYVV =    0.6294;
yVV =    0.6795;
omegamVV =    8.6453;
omegayVV =    6.2599;
oy_hf_initVV =    2.1589;
oy_hn_initVV =    1.2628;
oy_lf_initVV =    1.5749;
oy_ln_initVV =    1.4290;
om_hf_initVV =    2.9619;
om_hn_initVV =    1.3332;
om_lf_initVV =    2.8993;
om_ln_initVV =    1.4348;
N_forceVV =    1.4243;
margcostVV =   14.8932;
l_hftyVV =    0.2014;
l_hntyVV =    0.1259;
l_lftyVV =    0.1511;
l_lntyVV =    0.1385;
l_hftmVV =    0.2433;
l_hntmVV =    0.1145;
l_lftmVV =    0.2290;
l_lntmVV =    0.1217;
SVV =    0.0388;
mVV =    0.0553;
thigtnessVV =    2.0288;
VsVV =    0.0787;
%ejVV =    0.9178;
%y_piVV =   28.5354;