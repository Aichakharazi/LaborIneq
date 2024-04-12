function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = labormodel(isLog) 

if nargin == 0
    isLog = 0;
end



%Order of approximation desired 
approx = 2;

%Define parameters
syms de Sfe ej alp N_forcef e_m u_m ms thigtness Vs zz eeta bbeta aalpha 
syms ssp ssigmaparam xiparam ssigmaTFP rrhoTFP 
syms theta_hf theta_hn theta_lf theta_ln phi_hf phi_hn phi_lf phi_ln 
syms varrho upsilon varepsilon zetaparam ddeltaCov 


%Define variables 
syms zzV rsV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV 
syms LYV LMV l_hftmV l_hntmV l_lftmV l_lntmV l_hftyV l_hntyV l_lftyV l_lntyV 
syms yV margcostV lambdacV cmV cyV coldV e2V omegamV omegayV SV N_forceV e_mV u_mV mV thigtnessV VsV 
syms zzVV rsVV oy_hf_initVV oy_hn_initVV oy_lf_initVV oy_ln_initVV om_hf_initVV om_hn_initVV om_lf_initVV om_ln_initVV 
syms LYVV LMVV l_hftmVV l_hntmVV l_lftmVV l_lntmVV l_hftyVV l_hntyVV l_lftyVV l_lntyVV 
syms yVV margcostVV lambdacVV cmVV cyVV coldVV e2VV omegamVV omegayVV SVV N_forceVV e_mVV u_mVV mVV thigtnessVV VsVV 






% PRRODUCTIVITY SHOCK
%negative productivity shock 
%e27 = -(log(zzVV) - log(zz) ) + (1/zz)*  rrhoTFP * (log(zzV) - log(zz) )+  (1/zz)*(ddeltaCov * e2V) - (1/zz) * ssigmaTFP * e1; 

e1 = -(log(zzVV) - log(zz) ) + (1/zz)*  rrhoTFP * (log(zzV) - log(zz) )+  (1/zz)*(ddeltaCov * e2V); 
%e1 = -zzVV +   rrhoTFP * zzV  +  ddeltaCov * e2V; 

e2 = - e2VV  + e2V ;

% SV 
e3 = -SVV + SV;


% cy
e4 = cyV - 1 / ( lambdacVV * (1+rsV) ) ;
%e4 = cyhfV - 1 / ( lambdacVV * (1+rsV) ) ;
%e4 = cyhnV - 1 / ( lambdacVV * (1+rsV) ) ;
%e4 = cylfV - 1 / ( lambdacVV * (1+rsV) ) ;
%e4 = cylnV - 1 / ( lambdacVV * (1+rsV) ) ;


% r
e5 = lambdacV - lambdacVV*bbeta*(1+rsV) ; 

% cm
e6 = bbeta/cmV - lambdacV;
%e6 = bbeta/cmhfV - lambdacV;
%e6 = bbeta/cmhnV - lambdacV;
%e6 = bbeta/cmlfV - lambdacV;
%e6 = bbeta/cmlnV - lambdacV;

% co
e7 = bbeta^(2)/coldVV - lambdacV/(1+rsV); 

% LY
e8 = -LMV + ((lambdacV*e_mV*omegamV)/bbeta)^(1/eeta) ; 
% e8 = -l_hftmV + ((lambdacV*e_mV*om_hf_initV)/bbeta)^(1/eeta) ; 
% e8 = -l_hntmV + ((lambdacV*e_mV*om_hn_initV)/bbeta)^(1/eeta) ; 
% e8 = -l_lftmV + ((lambdacV*e_mV*om_lf_initV)/bbeta)^(1/eeta) ; 
% e8 = -l_lntmV + ((lambdacV*e_mV*om_ln_initV)/bbeta)^(1/eeta) ; 


% lm
e9 = -LYV + (lambdacV*e_mV*omegayV* (1+rsV) )^(1/eeta)   ; 
% e9 = -l_hftyV + (lambdacV*e_mV*oy_hf_initV* (1+rsV) )^(1/eeta)   ; 
% e9 = -l_hntyV + (lambdacV*e_mV*oy_hn_initV* (1+rsV) )^(1/eeta)   ; 
% e9 = -l_lftyV + (lambdacV*e_mV*oy_lf_initV* (1+rsV) )^(1/eeta)   ; 
% e9 = -l_lntyV + (lambdacV*e_mV*oy_ln_initV* (1+rsV) )^(1/eeta)   ; 





% euler equations and foc firm
%  omegayV  argcostV  omegamV
e10 = - margcostV +  (omegayV*((omegayV/omegamV)*((1-aalpha)/aalpha))^(ssp/(1-ssp)) + omegamV )/(zzV*(aalpha*((omegayV/omegamV*(1-aalpha)/aalpha)^(ssp/(1-ssp)))^(((ssp-1)/(ssp))) + (1-aalpha)* (1)^((ssp-1)/ssp)  )^(ssp/(1-ssp)));

%e10 = - omegayV /omegamV + aalpha/(1-aalpha)  * (LYV/LMV)^(1/ssp-1);
e11 = - omegayV + margcostV * aalpha *zzV * (LYV )^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);
e12 = - omegamV + margcostV * (1-aalpha) *zzV * (LMV )^(1/ssp-1)  * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);


% y
e13 = -yV + zzV * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1));
 



%    8 omegas  labor 
e14 = oy_hf_initV - margcostV * aalpha *zzV * ( l_hftyV)^(1/varrho-1) * phi_hf *  (( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV)^(1/varrho) * phi_hn )^(varrho/ssigmaparam-1) * ((( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lftyV)^(1/upsilon) * phi_lf + (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1)  *(LYV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e15 = oy_hn_initV - margcostV * aalpha *zzV * (l_hntyV)^(1/varrho-1) * phi_hn  * (( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV)^(1/varrho) * phi_hn )^(varrho/ssigmaparam-1)  * ((( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV   )^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lftyV )^(1/upsilon)  * phi_lf+ (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1) *(LYV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e16 = oy_lf_initV - margcostV * aalpha *zzV * (l_lftyV)^(1/upsilon-1) * phi_lf * ((l_lftyV)^(1/upsilon) * phi_lf + (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam-1) * ((( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lftyV)^(1/upsilon) * phi_lf + (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1) *(LYV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e17 = oy_ln_initV - margcostV * aalpha *zzV * (l_lntyV)^(1/upsilon-1) * phi_ln * ((l_lftyV)^(1/upsilon) * phi_lf + (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam-1) * ((( l_hftyV)^(1/varrho) * phi_hf + (l_hntyV)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lftyV)^(1/upsilon) * phi_lf + (l_lntyV)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1) *(LYV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e18 = om_hf_initV - margcostV * (1-aalpha) *zzV * ( l_hftmV)^(1/varepsilon-1) * theta_hf * (( l_hftmV)^(1/varepsilon) * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam-1) * ((( l_hftmV)^(1/varepsilon)  * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn )^(varepsilon/xiparam) + ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1) *(LMV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e19 = om_hn_initV -  margcostV * (1-aalpha) *zzV * (l_hntmV)^(1/varepsilon-1) * theta_hn * (( l_hftmV)^(1/varepsilon) * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam-1) * ((( l_hftmV)^(1/varepsilon) * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1) *(LMV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e20 = om_lf_initV -  margcostV * (1-aalpha) *zzV * (l_lftmV)^(1/zetaparam-1) * theta_lf * ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam-1) * ((( l_hftmV)^(1/varepsilon) * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1) *(LMV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);

e21 = om_ln_initV -  margcostV * (1-aalpha) *zzV * (l_lntmV)^(1/zetaparam-1) * theta_ln  * ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam-1) * ((( l_hftmV)^(1/varepsilon) * theta_hf + (l_hntmV)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftmV)^(1/zetaparam) * theta_lf + (l_lntmV)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1) *(LMV)^(1/ssp-1) * ( aalpha * (LYV)^(1/ssp) + (1-aalpha)* (LMV)^(1/ssp) )^(ssp-1);


% 8 labor  l
e22 = theta_hf * LYV -  l_hftyV; 
e23 = theta_hn * LYV -  l_hntyV;
e24 = theta_lf * LYV -  l_lftyV;
e25 = theta_ln * LYV -  l_lntyV;
e26 = phi_hf  * LMV  -  l_hftmV;
e27 = phi_hn  * LMV  -  l_hntmV;
e28 = phi_lf  * LMV  -  l_lftmV; 
e29 = phi_ln  * LMV  -  l_lntmV;


%  LYV =  l_hftyV +  l_hntyV + l_lftyV + l_lntyV;
%  LMV  =  l_hftmV +   l_hntmV + l_lftmV + l_lntmV;
 


% labor force
e30 = -N_forceVV + (1-de)* N_forceV  + ej *SV;

% employemnent 
e31 = -e_mV + (LYV+LMV)/N_forceV ;
% unemployment 
e32 = -u_mV + 1-e_mV;
% matching function
e33 = -mV + N_forceV *SV;
% labor market tightness
e34 = -thigtnessV + N_forceVV^(1/(1-alp));
% Vaccancies
e35 = -VsV + SV * thigtnessV;





e36 = -1/cyV + (LYV)^(eeta)/(e_mV*omegayV)  ;
% euler equations
%e10 = -1/cmV + (LMV)^(eeta)/(e_mV*omegamV);
%e11 = -1/cyV + (LYV)^(eeta)/(e_mV*omegayV);
%e12 = -bbeta/cmV + 1/(cyV*(1+rsV))  ;
% y_piV = bbeta*y_piV(+1)- kappa*margcostV; // philips curve





f = [e1;e2;e26;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25;e27;e28;e29;e30;e31;e32;e33;e34;e35;e36]
   
%f = [e27;e28]e8;e9;;e26
   
% Define the vector of controls, y, and states, x 


x  = [zzV N_forceV coldV ] 
        
y  = [lambdacV e2V SV cyV cmV  LYV LMV omegamV omegayV margcostV rsV yV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV l_hftmV l_hntmV l_lftmV l_lntmV l_hftyV l_hntyV l_lftyV l_lntyV  e_mV u_mV  mV thigtnessV VsV] 

xp = [zzVV N_forceVV coldVV] 

yp = [lambdacVV e2VV SVV cyVV cmVV  LYVV LMVV omegamVV omegayVV margcostVV rsVV yVV oy_hf_initVV oy_hn_initVV oy_lf_initVV oy_ln_initVV om_hf_initVV om_hn_initVV om_lf_initVV om_ln_initVV l_hftmVV l_hntmVV l_lftmVV l_lntmVV l_hftyVV l_hntyVV l_lftyVV l_lntyVV e_mVV u_mVV  mVV thigtnessVV VsVV] 
     
%yp = [rsVV oy_hf_initVV oy_hn_initVV oy_lf_initVV oy_ln_initVV om_hf_initVV om_hn_initVV om_lf_initVV om_ln_initVV LYVV LMVV l_hftmVV l_hntmVV l_lftmVV l_lntmVV l_hftyVV l_hntyVV l_lftyVV l_lntyVV yVV lambdacVV cmVV cyVV coldVV omegamVV omegayVV u_mVV e_mVV mVV N_forceVV VsVV thigtnessVV] 
       
       
%if isLog == 1
%    f = subs(f, [x, rsV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV LYV  LMV l_hftmV l_hntmV l_lftmV l_lntmV l_hftyV l_hntyV l_lftyV l_lntyV yV lambdacV cmV cyV coldV e2V omegamV omegayV u_mV e_mV mV N_forceV VsV thigtnessV, xp, rsVV oy_hf_initVV oy_hn_initVV oy_lf_initVV oy_ln_initVV om_hf_initVV om_hn_initVV om_lf_initVV om_ln_initVV LYVV LMVV l_hftmVV l_hntmVV l_lftmVV l_lntmVV l_hftyVV l_hntyVV l_lftyVV l_lntyVV yVV lambdacVV cmVV cyVV coldVV e2VV omegamVV omegayVV u_mVV e_mVV mVV N_forceVV VsVV thigtnessVV], ...
%        (exp([x, rsV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV LYV  LMV l_hftmV l_hntmV l_lftmV l_lntmV l_hftyV l_hntyV l_lftyV l_lntyV yV lambdacV cmV cyV coldV e2V omegamV omegayV u_mV e_mV mV N_forceV VsV thigtnessV, xp, rsVV oy_hf_initVV oy_hn_initVV oy_lf_initVV oy_ln_initVV om_hf_initVV om_hn_initVV om_lf_initVV om_ln_initVV LYVV LMVV l_hftmVV l_hntmVV l_lftmVV l_lntmVV l_hftyVV l_hntyVV l_lftyVV l_lntyVV yVV lambdacVV cmVV cyVV coldVV e2VV omegamVV omegayVV u_mVV e_mVV mVV N_forceVV VsVV thigtnessVV])));
%end


if isLog == 1
    f = subs(f, [x,y ,xp, yp ], ...
        (exp([x, y, xp,  yp])));
end


%if isLog == 1
%    f = subs(f, [x,  ], ...
%        (exp([x, ])));
%end

    


%Compute analytical derivatives of f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);
