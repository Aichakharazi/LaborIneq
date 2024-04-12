%% Model Equilibruim   Without capital 
%% Preferences 
eeta       = 1.5;          % 1.5the curvature on the disutility of labor  1-4
bbeta      = 0.95;       % disocunt factor 

%% Production function
aalpha      = 0.42;       % CPS youth labor share
% Aggregate shock z_{t} to be estimated
ssigmaTFP  = 0.08 ;        % standard deviation of TFP shock
rrhoTFP    = 0.98;       % quarterly autocorrelation of TFP shock
zz         = 1;

inzz = 1;

ddeltaCov = -0.0017;
%% Labor-labor elasticities : 

% elasticties (Euler equation estimation)
ssp_elas          =  0.076580  ;               % CPS estimated  (robust and with County and year FE)
ssp               = (1/(1+ssp_elas)) ;       % CPS
sigmaparam_elas   = 0.03924;                 % CPS estimated  (!!!! not robust and with County and year FE)
ssigmaparam       = (1/(sigmaparam_elas +1)); % CPS
xiparam_elas  = 0.07146;                       % CPS estimated  (robust and with County and year FE)
xiparam =(1/(xiparam_elas+1));               % CPS


% elasticties (Euler equation estimation)
varrho_elas      = 0.14092;                     % CPS estimated (robust  with and without County and year FE)
varrho           = 1/(varrho_elas +1);          % CPS 
upsilon_elas     = 0.11723 ;                    % CPS estimated (!!!!!not robust and with County and year FE)
upsilon          = 1/(upsilon_elas +1);         % CPS 
varepsilon_elas  = 0.05898;                      % CPS estimated (!!!!!not robust and with County and year FE)
varepsilon       = 1/(varepsilon_elas +1)  ;    % CPS 
zetaparam_elas   = 0.11216;                      % CPS estimated (not robust (significant only) without County and year FE)
zetaparam        = 1/(zetaparam_elas +1)  ;     % CPS

%% solve for interest rate and share of assets holding when young and middle aged
% interest rate
rs               = (1/bbeta)-1;
zma = 0.5;
zya = 0.5;

%% Household optimality conditions
% solve for cy cm co ly lm 
lambdac = 0.04012; % initial guess
cm = bbeta/lambdac;  
cy = 1/(lambdac*(1+rs));
co=(bbeta^(2)*(1+rs))/lambdac;

%% average wage rate
omegay=14.50;  % CPS (from the CPS year 2018)
omegam=13.36;  % CPS (from the CPS year 2018)
%% Household optimality conditions - labor
e_m = 0.95;   % employmnent rate

u_m = 1- e_m;  % unemployment rate

LM = ((lambdac *e_m * omegam)/bbeta)^(1/eeta); % lambdac = (bbeta*(LM)^(eeta))/omegam;
LY = (lambdac  *e_m * omegay*(1+rs) )^(1/eeta) ;





%% Labor : Distribution parameters  

% labor-type intensity is the ratio of labor j to \sigma j (total labor) in the production process
theta_hf  = 0.351;    % CPS high skill foreign labor share
theta_hn  = 0.256;    % CPS high skill native labor share
theta_lf  = 0.242;    % CPS low skill foreign labor share
theta_ln  = 0.249;    % CPS low skill native labor share


phi_hf    = 0.263;      % CPS high skill foreign labor share
phi_hn    = 0.265;      % CPS high skill native labor share
phi_lf    = 0.238;      % CPS low skill foreign labor share
phi_ln    = 0.232;      % CPS low skill native labor share

%% solve for labor
%LL =1;
%LY=aalpha*LL;
%LM = (1-aalpha)*LL;

l_hfty = phi_hf *LY;
l_hnty = phi_hn *LY;
l_lfty = phi_lf *LY;
l_lnty = phi_ln *LY;

l_hftm = theta_hf *LM;
l_hntm = theta_hn *LM;
l_lftm = theta_lf *LM;
l_lntm = theta_ln *LM;





%% solve for average wage
margcost = 1;

%%
  

% recommmneded 
%options = optimset('disp', 'iter', 'LargeScale', 'off', 'TolX', 1e-10, 'TolFun', 1e-12, 'MaxIter', 100000, 'MaxFunEvals', 1000000, 'Algorithm', 'trust-region-reflective');

options=optimset('disp','iter','LargeScale','off','TolX',1e-8,'TolFun',1e-10,'MaxIter',100000000,'MaxFunEvals',10000000000);

x0=[22.95; 16.15; 9.78; 6.47; 21.01; 16.26; 11.06; 8.45]; %

rex = lsqnonlin(@(x) solveomegas(x, LY, LM, aalpha,ssp, ssigmaparam, xiparam,zz, theta_hf, theta_hn, theta_lf, theta_ln, phi_hf,  phi_hn, phi_lf, phi_ln, varrho, upsilon, varepsilon, zetaparam,  l_hftm, l_hntm, l_lftm, l_lntm, l_hfty, l_hnty, l_lfty, l_lnty, margcost),x0,[],[],options);


% LY, LM, 

oy_hf_init= rex(1);
oy_hn_init = rex(2);
oy_lf_init = rex(3);
oy_ln_init= rex(4);
 
om_hf_init= rex(5);
om_hn_init= rex(6);
om_lf_init= rex(7);
om_ln_init= rex(8); 






%% Producer marginal cost:  
total_l_cost = l_hftm * om_hf_init  + l_hntm * om_hn_init + l_lftm * om_lf_init  + l_lntm * om_ln_init +l_hfty * oy_hf_init + l_hnty * oy_hn_init + l_lfty * oy_lf_init+ l_lnty * oy_ln_init  ;
y = zz * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp);
margcostgues = total_l_cost / y ;


%%
aa = y - cm - cy - co;



%% wage inflation and real wages


oy_hf_init_real = rex(1);
oy_hn_init_real  = rex(2);
oy_lf_init_real  = rex(3);
oy_ln_init_real = rex(4);
 
om_hf_init_real = rex(5);
om_hn_init_real = rex(6);
om_lf_init_real = rex(7);
om_ln_init_real = rex(8); 

% inflation assume it is positive  suppose it is equal to 0.001 instead of
% zero

oy_hf_init_pi  = oy_hf_init - oy_hf_init_real + 0.001;
oy_hn_init_pi  = oy_hn_init - oy_hn_init_real + 0.001;
oy_lf_init_pi  = oy_lf_init - oy_lf_init_real + 0.001;
oy_ln_init_pi  = oy_ln_init - oy_ln_init_real + 0.001;
 
om_hf_init_pi  = om_hf_init - om_hf_init_real + 0.001;
om_hn_init_pi  = om_hn_init - om_hn_init_real + 0.001;
om_lf_init_pi  = om_lf_init - om_lf_init_real + 0.001;
om_ln_init_pi  = om_ln_init - om_ln_init_real + 0.001; 




%%
Vrs = 0.06;
alp = 0.5;
%u_m = 0.05;
% Initial guess for parameters
x0 = [0.95,1.4,0.6,1.7];
options=optimset('disp','iter','LargeScale','off','TolX',1e-8,'TolFun',1e-10,'MaxIter',100000000,'MaxFunEvals',10000000000);

% Use lsqnonlin to find parameters that minimize the residuals
result = lsqnonlin(@(x) equationSystem(x,alp, Vrs,u_m, LY, LM),x0,[],[],options);
  
e_mf =result(1);
N_forcef = result(2);

Qsf = result(3);

tightnessf  = result(4);

% Display the result
disp('Optimal Parameters:');
disp(result);


%% save parameters and variables in steady states
save ss_model_LaborIneq eeta bbeta aalpha... 
ssp ssigmaparam xiparam...
ssigmaTFP rrhoTFP zz... 
theta_hf theta_hn theta_lf theta_ln...  
phi_hf  phi_hn phi_lf phi_ln...  
varrho upsilon varepsilon zetaparam...  
rs zma zya omegay omegam...
oy_hf_init oy_hn_init oy_lf_init oy_ln_init... 
om_hf_init om_hn_init om_lf_init om_ln_init... 
LY LM...
l_hftm l_hntm l_lftm l_lntm... 
l_hfty l_hnty l_lfty l_lnty... 
total_l_cost y  margcost... 
lambdac cm cy co aa margcostgues ddeltaCov inzz...
oy_hf_init_real oy_hn_init_real oy_lf_init_real oy_ln_init_real...
om_hf_init_real om_hn_init_real om_lf_init_real om_ln_init_real...
oy_hf_init_pi oy_hn_init_pi oy_lf_init_pi oy_ln_init_pi...
om_hf_init_pi om_hn_init_pi om_lf_init_pi om_ln_init_pi;




%% get some data from BEA-CPS (aggregated)









