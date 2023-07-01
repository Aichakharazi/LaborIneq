%% Model Equilibruim
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
ssp               = (-1/(ssp_elas)) ;       % CPS
sigmaparam_elas   = 0.03924;                 % CPS estimated  (!!!! not robust and with County and year FE)
ssigmaparam       = (1/(sigmaparam_elas+1)); % CPS
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

LM = ((lambdac*omegam)/bbeta)^(1/eeta); % lambdac = (bbeta*(LM)^(eeta))/omegam;
LY = (lambdac* omegay*(1+rs) )^(1/eeta) ;


%% Labor : Distribution parameters  

% labor-type intensity is the ratio of labor j to \sigma j (total labor) in the production process
% Middle-ages Labor supply

ttheta_hf  = 0.351;    % CPS high skill foreign labor share
ttheta_hn  = 0.256;    % CPS high skill native labor share
ttheta_lf  = 0.242;    % CPS low skill foreign labor share
ttheta_ln  = 0.249;    % CPS low skill native labor share

% Youth Labor supply

pphi_hf    = 0.263;      % CPS high skill foreign labor share
pphi_hn    = 0.265;      % CPS high skill native labor share
pphi_lf    = 0.238;      % CPS low skill foreign labor share
pphi_ln    = 0.232;      % CPS low skill native labor share

%% solve for labor
%LL =1;
%LY=aalpha*LL;
%LM = (1-aalpha)*LL;

LY_hf = pphi_hf *LY;
LY_hn = pphi_hn *LY;
LY_lf = pphi_lf *LY;
LY_ln = pphi_ln *LY;

LM_hf = ttheta_hf *LM;
LM_hn = ttheta_hn *LM;
LM_lf = ttheta_lf *LM;
LM_ln = ttheta_ln *LM;





%% solve for average wage
margcost = 1;

%%
  



options=optimset('disp','iter','LargeScale','off','TolX',1e-8,'TolFun',1e-10,'MaxIter',100000000,'MaxFunEvals',10000000000);

x0=[22.95; 16.15; 9.78; 6.47; 21.01; 16.26; 11.06; 8.45]; %

rex = lsqnonlin(@(x) solveomegas(x,aalpha,ssp, ssigmaparam, xiparam,zz, ttheta_hf, ttheta_hn, ttheta_lf, ttheta_ln, pphi_hf,  pphi_hn, pphi_lf, pphi_ln, varrho, upsilon, varepsilon, zetaparam, LY, LM, LM_hf, LM_hn, LM_lf, LM_ln, LY_hf, LY_hn, LY_lf, LY_ln, margcost),x0,[],[],options);



oy_hf_init= rex(1);
oy_hn_init = rex(2);
oy_lf_init = rex(3);
oy_ln_init= rex(4);
 
om_hf_init= rex(5);
om_hn_init= rex(6);
om_lf_init= rex(7);
om_ln_init= rex(8); 






%% Producer marginal cost:  
total_l_cost = LM_hf * om_hf_init  + LM_hn * om_hn_init + LM_lf * om_lf_init  + LM_ln * om_ln_init +LY_hf * oy_hf_init + LY_hn * oy_hn_init + LY_lf * oy_lf_init+ LY_ln * oy_ln_init  ;
y = zz * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1));
margcostgues = total_l_cost / y ;


%%
aa = y - cm - cy - co;



%%



%% save parameters and variables in steady states
save ss_model_LaborIneq eeta bbeta aalpha... 
ssp ssigmaparam xiparam...
ssigmaTFP rrhoTFP zz... 
ttheta_hf ttheta_hn ttheta_lf ttheta_ln...  
pphi_hf  pphi_hn pphi_lf pphi_ln...  
varrho upsilon varepsilon zetaparam...  
rs zma zya omegay omegam...
oy_hf_init oy_hn_init oy_lf_init oy_ln_init... 
om_hf_init om_hn_init om_lf_init om_ln_init... 
LY LM...
LM_hf LM_hn LM_lf LM_ln... 
LY_hf LY_hn LY_lf LY_ln... 
total_l_cost y  margcost... 
lambdac cm cy co aa margcostgues ddeltaCov inzz;




%% get some data from BEA-CPS (aggregated)






