// declare variables

var zzV, inzzV,
rsV, 
oy_hf_initV, oy_hn_initV, oy_lf_initV, oy_ln_initV, 
om_hf_initV, om_hn_initV, om_lf_initV, om_ln_initV,
LYV, LMV, 
LM_hfV, LM_hnV, LM_lfV, LM_lnV,
LY_hfV, LY_hnV, LY_lfV, LY_lnV,
yV, margcostV, 
lambdacV, cmV, cyV, coV,
yobs, 
labobs,
total_l_costV;


varexo e1, e2;

// declare parameteres (including variables in steady states)

parameters eeta, bbeta, aalpha, 
ssp, ssigmaparam, xiparam,
ssigmaTFP, rrhoTFP, zz, 
ttheta_hf, ttheta_hn, ttheta_lf, ttheta_ln,  
pphi_hf, pphi_hn, pphi_lf, pphi_ln, 
varrho, upsilon, varepsilon, zetaparam,  
rs, zma, zya, omegay, omegam,
oy_hf_init, oy_hn_init, oy_lf_init, oy_ln_init, 
om_hf_init, om_hn_init, om_lf_init, om_ln_init, 
LY, LM,
LM_hf, LM_hn, LM_lf, LM_ln, 
LY_hf, LY_hn, LY_lf, LY_ln,
total_l_cost, y, margcost, 
lambdac, cm, cy, co, aa, margcostgues, ddeltaCov, inzz;






load ss_model_LaborIneq.mat;
for i=1:length(M_.params)
    deep_parameter_name = M_.param_names(i,:);
    eval(['M_.params(i)  = ' deep_parameter_name ' ;'])
end;



model;

// Consumer Bloc
// F.O.C. wrt cy
//cyV - 1 / ( lambdacV(+1) * (1+rsV) ) ; // we shifted the timing of all variables
cyV(-1) - 1 / ( lambdacV * (1+rsV(-1)) ) ; //


// Consumer Euler equation 
lambdacV - lambdacV(+1)*bbeta*(1+rsV) ; 


// F.O.C. wrt cm
bbeta/cmV - lambdacV;

// F.O.C. wrt co
//bbeta^(2)/coV - lambdacV(-1)/(1+rsV(-1)); //  we shifted the timing of all variables
bbeta^(2)/coV(+1) - lambdacV/(1+rsV); // 


// F.O.C. wrt to lm 
-LMV + ((lambdacV*omegam)/bbeta)^(1/eeta); 


// F.O.C. wrt to ly
//-LYV + (lambdacV(+1)* omegay*(1+rsV) )^(1/eeta) ; //  we shifted the timing of all variables
//-LYV(-1) + (lambdacV* omegay(-1)*(1+rsV(-1)) )^(1/eeta) ; //  

-LYV + (lambdacV* omegay*(1+rsV(-1)) )^(1/eeta) ; //  



// producer Bloc

// production function
yV = zzV * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1));
 



total_l_costV = LM_hfV * om_hf_initV  + LM_hnV * om_hn_initV + LM_lfV * om_lf_initV  + LM_lnV * om_ln_initV +LY_hfV * oy_hf_initV + LY_hnV * oy_hn_initV + LY_lfV * oy_lf_initV+ LY_lnV * oy_ln_initV  ;
margcostgues = total_l_costV / yV ;



// F.O.C. wrt labor (same order of equations in the appendix) 
oy_hf_initV - margcostV * aalpha *zzV * ( LY_hfV * pphi_hf)^(1/varrho-1) *  (( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam-1) 
* ((( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1) 
*(LYV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

oy_hn_initV - margcostV * aalpha *zzV * (LY_hnV * pphi_hn )^(1/varrho-1) * (( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam-1) 
* ((( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)
*(LYV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


 oy_lf_initV - margcostV * aalpha *zzV * (LY_lfV * pphi_lf)^(1/upsilon-1) * ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam-1)
* ((( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)
*(LYV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


 oy_ln_initV - margcostV * aalpha *zzV * (LY_lnV * pphi_ln)^(1/upsilon-1) * ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam-1)
* ((( LY_hfV * pphi_hf)^(1/varrho) + (LY_hnV * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lfV * pphi_lf)^(1/upsilon) + (LY_lnV * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)
*(LYV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);



om_hf_initV - margcostV * (1-aalpha) *zzV * ( LM_hfV * ttheta_hf)^(1/varepsilon-1) * (( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam-1)
* ((( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)
*(LMV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

om_hn_initV -  margcostV * (1-aalpha) *zzV * (LM_hnV * ttheta_hn)^(1/varepsilon-1) * (( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam-1)
* ((( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)
*(LMV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


om_lf_initV -  margcostV * (1-aalpha) *zzV * (LM_lfV * ttheta_lf)^(1/zetaparam-1) * ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam-1)
* ((( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)
*(LMV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

om_ln_initV -  margcostV * (1-aalpha) *zzV * (LM_lnV * ttheta_ln)^(1/zetaparam-1)  * ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam-1)
* ((( LM_hfV * ttheta_hf)^(1/varepsilon) + (LM_hnV * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lfV * ttheta_lf)^(1/zetaparam) + (LM_lnV * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)
*(LMV)^((ssp-1)/ssp-1) * ( aalpha * (LYV)^((ssp-1)/ssp) + (1-aalpha)* (LMV)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);



// labor-labor intensity 
//ttheta_hf * LYV -  LY_hfV;
//ttheta_hn * LYV -  LY_hnV;
//ttheta_lf * LYV -  LY_lfV;
//ttheta_ln * LYV -  LY_lnV ;
//pphi_hf  * LMV  -  LM_hfV;
//pphi_hn  * LMV  -  LM_hnV;
//pphi_lf  * LMV  -  LM_lfV; 
//pphi_ln  * LMV  -  LM_lnV;


ttheta_hf * LYV -  LY_hfV;
ttheta_hn * LYV -  LY_hnV;
ttheta_lf * LYV -  LY_lfV;
ttheta_ln * LYV -  LY_lnV ;
pphi_hf  * LMV  -  LM_hfV;
pphi_hn  * LMV  -  LM_hnV;
pphi_lf  * LMV  -  LM_lfV; 
pphi_ln  * LMV  -  LM_lnV;



// Market clearing conditions 
aa = yV - cmV - cyV - coV;





// negative productivity shock 

//log(inzzV) - log(inzz) =  rrhoTFP * (log(inzzV(-1)) - log(inzz) )+  (1/inzz)*(ddeltaCov * e2) - (1/inzz) * ssigmaTFP * e1; 

//log(inzzV) - log(inzz) =  rrhoTFP * (log(inzzV(-1)) - log(inzz) ) - (1/inzz) * ssigmaTFP * e1; 
log(inzzV) - log(inzz) = (1/inzz)*  rrhoTFP * (log(inzzV(-1)) - log(inzz) )+  (1/inzz)*(ddeltaCov * e2) - (1/inzz) * ssigmaTFP * e1; 




-zzV +1/inzzV;








yobs      = log(yV) ;
labobs    = log(LMV + LYV) ;

end;




//var A;
//    periods 5, 6:9;
//    values 1.1, 1.05;
//



shocks ;
var e1; stderr 1 ; // productivity shock (to be estimated) 
var e2; stderr 1; // affected share shock (to be estimated) 
end;





initval;
inzzV = inzz;
zzV = zz;
rsV =  rs; 
oy_hf_initV = oy_hf_init;
oy_hn_initV = oy_hn_init;  
oy_lf_initV = oy_lf_init;
oy_ln_initV = oy_ln_init;
om_hf_initV = om_hf_init;
om_hn_initV = om_hn_init;
om_lf_initV = om_lf_init;
om_ln_initV = om_ln_init;
LYV = LY;
LMV = LM; 
LM_hfV = LM_hf;
LM_hnV =  LM_hn; 
LM_lfV = LM_lf;
LM_lnV = LM_ln;
LY_hfV = LY_hf;
LY_hnV = LY_hn;
LY_lfV =  LY_lf; 
LY_lnV = LY_ln;
yV = y;
margcostV = margcost; 
lambdacV = lambdac; 
cmV = cm;
cyV = cy;
coV = co;
total_l_costV = total_l_cost;

//yobs      = log(y) ;
//labobs    = log(LM + LY) ;
end;



resid(1);
steady(solve_algo=2,maxit=1000);


check;

steady;
resid(1);



varobs yobs labobs;
//calib_smoother(datafile=dataqdynarefin,first_obs=11)yobs labobs;
//calib_smoother(datafile=dataqdynarefin,first_obs=11, filtered_vars, filter_step_ahead = [3:4])yobs;// labobs;



calib_smoother(datafile=macrotimeseries1, filtered_vars, filter_step_ahead = [3:4]) zzV rsV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV LYV LMV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV yV margcostV lambdacV cmV cyV coV total_l_costV;

//calib_smoother(datafile=dataqdynarefin, filtered_vars, filter_step_ahead = [3:4]) zzV rsV oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV om_hf_initV om_hn_initV om_lf_initV om_ln_initV LYV LMV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV yV margcostV lambdacV cmV cyV coV total_l_costV;

// macrotimeseries1


// Impulse responses functions ,simul_replic=100000
stoch_simul(order=1,pruning,k_order_solver,irf=80,simul_replic=10000000);

// Save results in file LaborIneq_model
LaborIneq_model=oo_.irfs;
save LaborIneq_model







T=80;

// plot IRF
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,LM_hfV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LY_hfV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);

//title('Labor - Foreign - High skill','FontSize',20)
t = ' A. Labor - Foreign - High skill ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,2)
plot(1:T,LM_hnV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LY_hnV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Labor - Native - High skill','FontSize',20)
t = ' B. Labor - Native - High skill  ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,3)
plot(1:T,LM_lfV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LY_lfV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Labor - Foreign - Low skill','FontSize',20)
t = ' C. Labor - Foreign - Low skill  ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,4)
plot(1:T,LM_lnV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on
plot(1:T,LY_lnV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Labor - Native - Low skill ','FontSize',20)
t = ' D. Labor - Native - Low skill  ';
title(t,'interpreter','latex')


axis tight; 
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

// the shock: affect negatively wages and hours increases



// plot IRF for wages
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,om_hf_initV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,oy_hf_initV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Wages - Foreign - High skill','FontSize',20)
t = ' A. Wages - Foreign - High skill  ';
title(t,'interpreter','latex')



axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,2)
plot(1:T,om_hn_initV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,oy_hn_initV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Wages - Native - High skill','FontSize',20)
t = ' B. Wages - Native - High skill  ';
title(t,'interpreter','latex')


axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,3)
plot(1:T,om_lf_initV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,oy_lf_initV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Wages - Foreign - Low skill','FontSize',20)
t = ' C. Wages - Foreign - Low skill  ';
title(t,'interpreter','latex')


axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,4)
plot(1:T,om_ln_initV_e1(1:T)*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on
plot(1:T,oy_ln_initV_e1(1:T)*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
//title('Wages - Native - Low skill ','FontSize',20)
t = ' D. Wages - Native - Low skill  ';
title(t,'interpreter','latex')

axis tight; 
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'Middle aged','Young'},'FontSize',20,'Location','best')


// Second Results about sensitivity of wages to covid shock (as defined in our model) 
// Benchmark model ddeltaCov = -0.0017 (responses of productivity to covid)
set_param_value('ddeltaCov', -0.0017);  // you don't nedd to specify it is already declared
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV; 

oy_hf_initV_var_low_response = oo_.var(strmatch('oy_hf_initV',var_list_,'exact'),strmatch('oy_hf_initV',var_list_,'exact'));
oy_hn_initV_var_low_response = oo_.var(strmatch('oy_hn_initV',var_list_,'exact'),strmatch('oy_hn_initV',var_list_,'exact'));
oy_lf_initV_var_low_response = oo_.var(strmatch('oy_lf_initV',var_list_,'exact'),strmatch('oy_lf_initV',var_list_,'exact'));
oy_ln_initV_var_low_response = oo_.var(strmatch('oy_ln_initV',var_list_,'exact'),strmatch('oy_ln_initV',var_list_,'exact'));

om_hf_initV_var_low_response = oo_.var(strmatch('om_hf_initV',var_list_,'exact'),strmatch('om_hf_initV',var_list_,'exact'));
om_hn_initV_var_low_response = oo_.var(strmatch('om_hn_initV',var_list_,'exact'),strmatch('om_hn_initV',var_list_,'exact'));
om_lf_initV_var_low_response = oo_.var(strmatch('om_lf_initV',var_list_,'exact'),strmatch('om_lf_initV',var_list_,'exact'));
om_ln_initV_var_low_response = oo_.var(strmatch('om_ln_initV',var_list_,'exact'),strmatch('om_ln_initV',var_list_,'exact'));




LY_hfV_var_low_response = oo_.var(strmatch('LY_hfV',var_list_,'exact'),strmatch('LY_hfV',var_list_,'exact'));
LY_hnV_var_low_response = oo_.var(strmatch('LY_hnV',var_list_,'exact'),strmatch('LY_hnV',var_list_,'exact'));
LY_lfV_var_low_response = oo_.var(strmatch('LY_lfV',var_list_,'exact'),strmatch('LY_lfV',var_list_,'exact'));
LY_lnV_var_low_response = oo_.var(strmatch('LY_lnV',var_list_,'exact'),strmatch('LY_lnV',var_list_,'exact'));

LM_hfV_var_low_response = oo_.var(strmatch('LM_hfV',var_list_,'exact'),strmatch('LM_hfV',var_list_,'exact'));
LM_hnV_var_low_response = oo_.var(strmatch('LM_hnV',var_list_,'exact'),strmatch('LM_hnV',var_list_,'exact'));
LM_lfV_var_low_response = oo_.var(strmatch('LM_lfV',var_list_,'exact'),strmatch('LM_lfV',var_list_,'exact'));
LM_lnV_var_low_response = oo_.var(strmatch('LM_lnV',var_list_,'exact'),strmatch('LM_lnV',var_list_,'exact'));



// if we increase the responses of productivity to the covid shock how wages will react
set_param_value('ddeltaCov',-0.2);
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;

oy_hf_initV_var_high_response = oo_.var(strmatch('oy_hf_initV',var_list_,'exact'),strmatch('oy_hf_initV',var_list_,'exact'));
oy_hn_initV_var_high_response = oo_.var(strmatch('oy_hn_initV',var_list_,'exact'),strmatch('oy_hn_initV',var_list_,'exact'));
oy_lf_initV_var_high_response = oo_.var(strmatch('oy_lf_initV',var_list_,'exact'),strmatch('oy_lf_initV',var_list_,'exact'));
oy_ln_initV_var_high_response = oo_.var(strmatch('oy_ln_initV',var_list_,'exact'),strmatch('oy_ln_initV',var_list_,'exact'));
om_hf_initV_var_high_response = oo_.var(strmatch('om_hf_initV',var_list_,'exact'),strmatch('om_hf_initV',var_list_,'exact'));
om_hn_initV_var_high_response = oo_.var(strmatch('om_hn_initV',var_list_,'exact'),strmatch('om_hn_initV',var_list_,'exact'));
om_lf_initV_var_high_response = oo_.var(strmatch('om_lf_initV',var_list_,'exact'),strmatch('om_lf_initV',var_list_,'exact'));
om_ln_initV_var_high_response = oo_.var(strmatch('om_ln_initV',var_list_,'exact'),strmatch('om_ln_initV',var_list_,'exact'));



LY_hfV_var_high_response = oo_.var(strmatch('LY_hfV',var_list_,'exact'),strmatch('LY_hfV',var_list_,'exact'));
LY_hnV_var_high_response = oo_.var(strmatch('LY_hnV',var_list_,'exact'),strmatch('LY_hnV',var_list_,'exact'));
LY_lfV_var_high_response = oo_.var(strmatch('LY_lfV',var_list_,'exact'),strmatch('LY_lfV',var_list_,'exact'));
LY_lnV_var_high_response = oo_.var(strmatch('LY_lnV',var_list_,'exact'),strmatch('LY_lnV',var_list_,'exact'));
LM_hfV_var_high_response = oo_.var(strmatch('LM_hfV',var_list_,'exact'),strmatch('LM_hfV',var_list_,'exact'));
LM_hnV_var_high_response = oo_.var(strmatch('LM_hnV',var_list_,'exact'),strmatch('LM_hnV',var_list_,'exact'));
LM_lfV_var_high_response = oo_.var(strmatch('LM_lfV',var_list_,'exact'),strmatch('LM_lfV',var_list_,'exact'));
LM_lnV_var_high_response = oo_.var(strmatch('LM_lnV',var_list_,'exact'),strmatch('LM_lnV',var_list_,'exact'));


set_param_value('ddeltaCov',0);
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;


oy_hf_initV_var_zero_response = oo_.var(strmatch('oy_hf_initV',var_list_,'exact'),strmatch('oy_hf_initV',var_list_,'exact'));
oy_hn_initV_var_zero_response = oo_.var(strmatch('oy_hn_initV',var_list_,'exact'),strmatch('oy_hn_initV',var_list_,'exact'));
oy_lf_initV_var_zero_response = oo_.var(strmatch('oy_lf_initV',var_list_,'exact'),strmatch('oy_lf_initV',var_list_,'exact'));
oy_ln_initV_var_zero_response = oo_.var(strmatch('oy_ln_initV',var_list_,'exact'),strmatch('oy_ln_initV',var_list_,'exact'));
om_hf_initV_var_zero_response = oo_.var(strmatch('om_hf_initV',var_list_,'exact'),strmatch('om_hf_initV',var_list_,'exact'));
om_hn_initV_var_zero_response = oo_.var(strmatch('om_hn_initV',var_list_,'exact'),strmatch('om_hn_initV',var_list_,'exact'));
om_lf_initV_var_zero_response = oo_.var(strmatch('om_lf_initV',var_list_,'exact'),strmatch('om_lf_initV',var_list_,'exact'));
om_ln_initV_var_zero_response = oo_.var(strmatch('om_ln_initV',var_list_,'exact'),strmatch('om_ln_initV',var_list_,'exact'));



LY_hfV_var_zero_response = oo_.var(strmatch('LY_hfV',var_list_,'exact'),strmatch('LY_hfV',var_list_,'exact'));
LY_hnV_var_zero_response = oo_.var(strmatch('LY_hnV',var_list_,'exact'),strmatch('LY_hnV',var_list_,'exact'));
LY_lfV_var_zero_response = oo_.var(strmatch('LY_lfV',var_list_,'exact'),strmatch('LY_lfV',var_list_,'exact'));
LY_lnV_var_zero_response = oo_.var(strmatch('LY_lnV',var_list_,'exact'),strmatch('LY_lnV',var_list_,'exact'));
LM_hfV_var_zero_response = oo_.var(strmatch('LM_hfV',var_list_,'exact'),strmatch('LM_hfV',var_list_,'exact'));
LM_hnV_var_zero_response = oo_.var(strmatch('LM_hnV',var_list_,'exact'),strmatch('LM_hnV',var_list_,'exact'));
LM_lfV_var_zero_response = oo_.var(strmatch('LM_lfV',var_list_,'exact'),strmatch('LM_lfV',var_list_,'exact'));
LM_lnV_var_zero_response = oo_.var(strmatch('LM_lnV',var_list_,'exact'),strmatch('LM_lnV',var_list_,'exact'));



// create a table for the results (when we increase the responses of productivity to covid shock , wages become more senstive and differ by type of worker)
fprintf('Wages omega_{hf}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(oy_hf_initV_var_low_response),sqrt(oy_hf_initV_var_high_response),sqrt(oy_hf_initV_var_zero_response));
fprintf('Wages omega_{hn}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(oy_hn_initV_var_low_response),sqrt(oy_hn_initV_var_high_response),sqrt(oy_hn_initV_var_zero_response));
fprintf('Wages omega_{lf}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(oy_lf_initV_var_low_response),sqrt(oy_lf_initV_var_high_response),sqrt(oy_lf_initV_var_zero_response));
fprintf('Wages omega_{ln}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(oy_ln_initV_var_low_response),sqrt(oy_ln_initV_var_high_response),sqrt(oy_ln_initV_var_zero_response));
fprintf('Wages omega_{hf}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(om_hf_initV_var_low_response),sqrt(om_hf_initV_var_high_response),sqrt(om_hf_initV_var_zero_response));
fprintf('Wages omega_{hn}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(om_hn_initV_var_low_response),sqrt(om_hn_initV_var_high_response),sqrt(om_hn_initV_var_zero_response));
fprintf('Wages omega_{lf}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(om_lf_initV_var_low_response),sqrt(om_lf_initV_var_high_response),sqrt(om_lf_initV_var_zero_response));
fprintf('Wages omega_{ln}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(om_ln_initV_var_low_response),sqrt(om_ln_initV_var_high_response),sqrt(om_ln_initV_var_zero_response));



fprintf('Hours L_{hf}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LY_hfV_var_low_response),sqrt(LY_hfV_var_high_response),sqrt(LY_hfV_var_zero_response));
fprintf('Hours L_{hn}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LY_hnV_var_low_response),sqrt(LY_hnV_var_high_response),sqrt(LY_hnV_var_zero_response));
fprintf('Hours L_{lf}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LY_lfV_var_low_response),sqrt(LY_lfV_var_high_response),sqrt(LY_lfV_var_zero_response));
fprintf('Hours L_{ln}^{y} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LY_lnV_var_low_response),sqrt(LY_lnV_var_high_response),sqrt(LY_lnV_var_zero_response));
fprintf('Hours L_{hf}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LM_hfV_var_low_response),sqrt(LM_hfV_var_high_response),sqrt(LM_hfV_var_zero_response));
fprintf('Hours L_{hn}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LM_hnV_var_low_response),sqrt(LM_hnV_var_high_response),sqrt(LM_hnV_var_zero_response));
fprintf('Hours L_{lf}^{m} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LM_lfV_var_low_response),sqrt(LM_lfV_var_high_response),sqrt(LM_lfV_var_zero_response));
fprintf('Hours L_{ln}^{n} Standard Deviation: \t %4.3f \t %4.3f \t %4.3f\n',sqrt(LM_lnV_var_low_response),sqrt(LM_lnV_var_high_response),sqrt(LM_lnV_var_zero_response));



// wage differential 
// reset the value of ddeltaCov  (as in benchmark model)
// now compute the wage differential 
set_param_value('ddeltaCov', -0.0017);
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;

// Foreign vs Native (with same caharcteristic)
wage_diff_oy_hf_initV_oy_hn_initV  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;


fprintf('Wages differential (omega_{hf}^{y} - omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV    );
fprintf('Wages differential (omega_{lf}^{y} - omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV    );
fprintf('Wages differential (omega_{hf}^{m} - omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV    );
fprintf('Wages differential (omega_{lf}^{m} - omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV    );



// young vs milddle-aged (with same caharcteristic)
wage_diff_oy_hf_initV_om_hf_initV  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_om_lf_initV  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_hn_initV_om_hn_initV  = (oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_ln_initV_om_ln_initV  = (oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;


fprintf('Wages differential (omega_{hf}^{y} - omega_{hf}^{m}) : \t %4.3f\n', wage_diff_oy_hf_initV_om_hf_initV    );
fprintf('Wages differential (omega_{lf}^{y} - omega_{lf}^{m}) : \t %4.3f\n', wage_diff_oy_lf_initV_om_lf_initV    );
fprintf('Wages differential (omega_{hn}^{y} - omega_{hn}^{m}) : \t %4.3f\n', wage_diff_oy_hn_initV_om_hn_initV    );
fprintf('Wages differential (omega_{ln}^{y} - omega_{ln}^{m}) : \t %4.3f\n', wage_diff_oy_ln_initV_om_ln_initV    );




//  high vs low skill  (with same caharcteristic)
wage_diff_oy_hf_initV_oy_lf_initV  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_lf_initV  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_hn_initV_oy_ln_initV  = (oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hn_initV_om_ln_initV  = (oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;


fprintf('Wages differential (omega_{hf}^{y} - omega_{lf}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_lf_initV    );
fprintf('Wages differential (omega_{hf}^{m} - omega_{lf}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_lf_initV    );
fprintf('Wages differential (omega_{hn}^{y} - omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_hn_initV_oy_ln_initV    );
fprintf('Wages differential (omega_{hn}^{m} - omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_hn_initV_om_ln_initV    );



// Wages premium


wage_diff_oy_hf_initV_oy_hn_initV2  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV2  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV2  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV2  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;


fprintf('Wages premium (omega_{hf}^{y} / omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV2    );
fprintf('Wages premium (omega_{lf}^{y} / omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV2    );
fprintf('Wages premium (omega_{hf}^{m} / omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV2    );
fprintf('Wages premium (omega_{lf}^{m} / omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV2    );






set_param_value('ddeltaCov',-0.2);
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;

wage_diff_oy_hf_initV_oy_hn_initV_ScenarioB  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV_ScenarioB  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV_ScenarioB  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV_ScenarioB  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;



fprintf('Wages differential (omega_{hf}^{y} - omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV_ScenarioB    );
fprintf('Wages differential (omega_{lf}^{y} - omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV_ScenarioB    );
fprintf('Wages differential (omega_{hf}^{m} - omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV_ScenarioB    );
fprintf('Wages differential (omega_{lf}^{m} - omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV_ScenarioB    );

wage_diff_oy_hf_initV_oy_hn_initV_ScenarioB2  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV_ScenarioB2  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV_ScenarioB2  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV_ScenarioB2  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;

fprintf('Wages premium (omega_{hf}^{y} / omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV_ScenarioB2    );
fprintf('Wages premium (omega_{lf}^{y} / omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV_ScenarioB2    );
fprintf('Wages premium (omega_{hf}^{m} / omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV_ScenarioB2    );
fprintf('Wages premium (omega_{lf}^{m} / omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV_ScenarioB2    );



//
//
//
//

set_param_value('ddeltaCov',0);
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;

wage_diff_oy_hf_initV_oy_hn_initV_ScenarioC  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV_ScenarioC  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV_ScenarioC  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV_ScenarioC  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) - oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;



fprintf('Wages differential (omega_{hf}^{y} - omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV_ScenarioC    );
fprintf('Wages differential (omega_{lf}^{y} - omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV_ScenarioC    );
fprintf('Wages differential (omega_{hf}^{m} - omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV_ScenarioC    );
fprintf('Wages differential (omega_{lf}^{m} - omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV_ScenarioC    );

wage_diff_oy_hf_initV_oy_hn_initV_ScenarioC2  = (oo_.steady_state(strmatch('oy_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_oy_lf_initV_oy_ln_initV_ScenarioC2  = (oo_.steady_state(strmatch('oy_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('oy_ln_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_hf_initV_om_hn_initV_ScenarioC2  = (oo_.steady_state(strmatch('om_hf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_hn_initV',M_.endo_names,'exact')) ) ;
wage_diff_om_lf_initV_om_ln_initV_ScenarioC2  = (oo_.steady_state(strmatch('om_lf_initV',M_.endo_names,'exact')) / oo_.steady_state(strmatch('om_ln_initV',M_.endo_names,'exact')) ) ;

fprintf('Wages premium (omega_{hf}^{y} / omega_{hn}^{y}) : \t %4.3f\n', wage_diff_oy_hf_initV_oy_hn_initV_ScenarioC2    );
fprintf('Wages premium (omega_{lf}^{y} / omega_{ln}^{y}) : \t %4.3f\n', wage_diff_oy_lf_initV_oy_ln_initV_ScenarioC2    );
fprintf('Wages premium (omega_{hf}^{m} / omega_{hn}^{m}) : \t %4.3f\n', wage_diff_om_hf_initV_om_hn_initV_ScenarioC2    );
fprintf('Wages premium (omega_{lf}^{m} / omega_{ln}^{m}) : \t %4.3f\n', wage_diff_om_lf_initV_om_ln_initV_ScenarioC2    );




//
//
//


// wage differential is not sensitive to the higher responses of productivity to covid shock
// We reset again the parameter value of the benchmark model
set_param_value('ddeltaCov',-0.0017); // reset to banchmark
stoch_simul(irf=0,order=1,nograph,noprint) oy_hf_initV oy_hn_initV oy_lf_initV oy_ln_initV  om_hf_initV om_hn_initV om_lf_initV om_ln_initV LM_hfV LM_hnV LM_lfV LM_lnV LY_hfV LY_hnV LY_lfV LY_lnV;

// we can do also 
// Young vs Middle aged


// wage differentail 





// plot IRF for wages
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,1,1)
plot(1:T,(om_hn_initV_e1(1:T)-om_hf_initV_e1(1:T))*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,(om_ln_initV_e1(1:T)-om_lf_initV_e1(1:T))*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Foreign - High skill','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'mhn  mhf','mln mlf'},'FontSize',20,'Location','best')

subplot(2,1,2)
plot(1:T,(oy_hn_initV_e1(1:T)-oy_hf_initV_e1(1:T))*100,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,(oy_ln_initV_e1(1:T)-oy_lf_initV_e1(1:T))*100,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Native - High skill','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'yhn yhf','yln ylf'},'FontSize',20,'Location','best')



//

set_param_value('ddeltaCov',-0.0017);
stoch_simul(order=1,irf=11,nograph,nomoments,nocorr,nofunctions);

// simulation scenarios 

// Exercice III : COUNTERFACTUALS  plots

innovations = oo_.SmoothedShocks(:,1);
n_points=232;

// first  


//set_param_value('ddeltaCov',0); //
//stoch_simul(order=1,irf=11,nograph,nomoments,nocorr,nofunctions);



// initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(n_points,M_.exo_nbr);      %create shock matrix with number of time periods in colums

// set all shocks for e1 and e2;
shock_matrix(:,strmatch('e1',M_.exo_names,'exact')) = innovations.e1';
shock_matrix(:,strmatch('e2',M_.exo_names,'exact')) = innovations.e2';
 

y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_with_covidshock = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     % deviation from steady state  


//set_param_value('ddeltaCov',0); //
//stoch_simul(order=1,irf=11,nograph,nomoments,nocorr,nofunctions);
        
// set productivity shocks only e1 ;
shock_matrix(:,strmatch('e1',M_.exo_names,'exact')) = innovations.e1';
shock_matrix(:,strmatch('e2',M_.exo_names,'exact')) = zeros(1,n_points);
 

y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_without_covidshock = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);         
        

LM_hfV_exp1 = y_with_covidshock(14,:);
LM_hfV_exp2 = y_without_covidshock(14,:);

LM_hnV_exp1 = y_with_covidshock(15,:);
LM_hnV_exp2 = y_without_covidshock(15,:);

LM_lfV_exp1 = y_with_covidshock(16,:);
LM_lfV_exp2 = y_without_covidshock(16,:);

LM_lnV_exp1 = y_with_covidshock(17,:);
LM_lnV_exp2 = y_without_covidshock(17,:);

LY_hfV_exp1 = y_with_covidshock(18,:);
LY_hfV_exp2 = y_without_covidshock(18,:);


LY_hnV_exp1 = y_with_covidshock(19,:);
LY_hnV_exp2 = y_without_covidshock(19,:);

LY_lfV_exp1 = y_with_covidshock(20,:);
LY_lfV_exp2 = y_without_covidshock(20,:);

LY_lnV_exp1 = y_with_covidshock(21,:);
LY_lnV_exp2 = y_without_covidshock(21,:);



y_V_exp1 = y_with_covidshock(22,:);
y_V_exp2 = y_without_covidshock(22,:);




T=232;
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,LM_hfV_exp1,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LM_hfV_exp2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Foreign - High skill - Young','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')


subplot(2,2,2)
plot(1:T,LM_hnV_exp1,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LM_hnV_exp2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Foreign - High skill - Young','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')


subplot(2,2,3)
plot(1:T,LM_lfV_exp1,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LM_lfV_exp2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Foreign - High skill - Young','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')

subplot(2,2,4)
plot(1:T,LM_lnV_exp1,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,LM_lnV_exp2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wages - Foreign - High skill - Young','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')









// 


T=232;
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,y_V_exp1,'LineWidth',3,'Color', [242, 109, 116]/255); hold on;
plot(1:T,y_V_exp2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Output','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')
