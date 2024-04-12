

isLog = 1;


%% 1) Model Solution
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f]= labormodel(isLog);





%Assign values to parameters and steady-state variables

[de, Sfe, ej, alp, N_forcef, e_m, u_m, ms, thigtness, Vs, zz, eeta, bbeta, aalpha, ssp, ssigmaparam, xiparam, ssigmaTFP, rrhoTFP, theta_hf, theta_hn, theta_lf, theta_ln, phi_hf, phi_hn, phi_lf, phi_ln, varrho, upsilon, varepsilon, zetaparam, ddeltaCov, zzV, rsV, oy_hf_initV, oy_hn_initV, oy_lf_initV, oy_ln_initV, om_hf_initV, om_hn_initV, om_lf_initV, om_ln_initV, LYV, LMV, l_hftmV, l_hntmV, l_lftmV, l_lntmV, l_hftyV, l_hntyV, l_lftyV, l_lntyV, yV, margcostV, lambdacV, cmV, cyV, coldV, e2V, omegamV, omegayV, SV, N_forceV, e_mV, u_mV, mV, thigtnessV, VsV, zzVV, rsVV, oy_hf_initVV, oy_hn_initVV, oy_lf_initVV, oy_ln_initVV, om_hf_initVV, om_hn_initVV, om_lf_initVV, om_ln_initVV, LYVV, LMVV, l_hftmVV, l_hntmVV, l_lftmVV, l_lntmVV, l_hftyVV, l_hntyVV, l_lftyVV, l_lntyVV, yVV, margcostVV, lambdacVV, cmVV, cyVV, coldVV, e2VV, omegamVV, omegayVV, SVV, N_forceVV, e_mVV, u_mVV, mVV, thigtnessVV, VsVV] = sstatelabor;
approx = 2;
eta = [1 0 0]'; % ssigmaTFP  1 1
%Compute numerical derivatives of f
num_eval


%% 2) First Order Approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp)



%% 4) Impulse response function




x0=zeros(size(hx,1),1);
x0(1)=-1; % negative shock to productivity

%%
T=40;
IRF = irf_1st(gx, hx, eta,T,x0);


%%





l_hftmV_e1  =IRF(:,20);
l_hntmV_e1  =IRF(:,21);
l_lftmV_e1  =IRF(:,22);
l_lntmV_e1  =IRF(:,23);

l_hftyV_e1  =IRF(:,24);
l_hntyV_e1  =IRF(:,25);
l_lftyV_e1  =IRF(:,26);
l_lntyV_e1  =IRF(:,27);


oy_hf_initV_e1    =IRF(:,12);
oy_hn_initV_e1    =IRF(:,13);
oy_lf_initV_e1    =IRF(:,14);
oy_ln_initV_e1    =IRF(:,15);

om_hf_initV_e1  =IRF(:,16);
om_hn_initV_e1  =IRF(:,17);
om_lf_initV_e1  =IRF(:,18);
om_ln_initV_e1  =IRF(:,19);

%%  plot IRF labor demand



T=20;


figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,l_hftmV_e1(1:T)*10000,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,l_hftyV_e1(1:T)*10000,'LineWidth',3,'Color', '#0F5C8C'); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);

t = ' A. Labor - Foreign - High skill ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-0.1 1.3])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')





subplot(2,2,2)
plot(1:T,l_hntmV_e1(1:T)*10000,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,l_hntyV_e1(1:T)*10000,'LineWidth',3,'Color', '#0F5C8C'); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Labor - Native - High skill','FontSize',20)
t = ' B. Labor - Native - High skill  ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-0.1 1.3])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,3)
plot(1:T,l_lftmV_e1(1:T)*10000,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,l_lftyV_e1(1:T)*10000,'LineWidth',3,'Color', '#0F5C8C'); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Labor - Foreign - Low skill','FontSize',20)
t = ' C. Labor - Foreign - Low skill  ';
title(t,'interpreter','latex')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-0.1 1.3])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,4)
plot(1:T,l_lntmV_e1(1:T)*10000,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on
plot(1:T,l_lntyV_e1(1:T)*10000,'LineWidth',3,'Color', '#0F5C8C'); hold on
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Labor - Native - Low skill ','FontSize',20)
t = ' D. Labor - Native - Low skill  ';
title(t,'interpreter','latex')


axis tight; 
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-0.1 1.3])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

% the shock: affect negatively wages and hours increases


export_fig irfhours.jpg  -m2.5 -painters


%% plot IRF for wages
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,2,1)
plot(1:T,om_hf_initV_e1(1:T)*100,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,oy_hf_initV_e1(1:T)*100,'LineWidth',3,'Color', '#0F5C8C'); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Wages - Foreign - High skill','FontSize',20)
t = ' A. Wages - Foreign - High skill  ';
title(t,'interpreter','latex')



axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-5 0.5])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')


subplot(2,2,2)
plot(1:T,om_hn_initV_e1(1:T)*100,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,oy_hn_initV_e1(1:T)*100,'LineWidth',3,'Color', '#0F5C8C'); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Wages - Native - High skill','FontSize',20)
t = ' B. Wages - Native - High skill  ';
title(t,'interpreter','latex')


axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-5 0.5])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,3)
plot(1:T,om_lf_initV_e1(1:T)*100,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,oy_lf_initV_e1(1:T)*100,'LineWidth',3,'Color', '#0F5C8C'); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Wages - Foreign - Low skill','FontSize',20)
t = ' C. Wages - Foreign - Low skill  ';
title(t,'interpreter','latex')


axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-5 0.5])
legend({'Middle aged','Young'},'FontSize',20,'Location','best')

subplot(2,2,4)
plot(1:T,om_ln_initV_e1(1:T)*100,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on
plot(1:T,oy_ln_initV_e1(1:T)*100,'LineWidth',3,'Color', '#0F5C8C'); hold on
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
%title('Wages - Native - Low skill ','FontSize',20)
t = ' D. Wages - Native - Low skill  ';
title(t,'interpreter','latex')

axis tight; 
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
%ylim([-5 0.5])

legend({'Middle aged','Young'},'FontSize',20,'Location','best')

export_fig irfwages.jpg  -m2.5 -painters
%%





%% Model Simulation (High kappa)

sig = 0.01;
num_sim = 10000;
e = randn(num_sim, 1);

x0 = [0];
[Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);

smpl = 0.9*num_sim:num_sim;


%% this is important
%%
%% Skip this (change kappa and run the code and then run this part)




x0=zeros(size(hx,1),1);
x0(1)=1; % shock to productivity



T=40;
IRFLow = irf_1st(gx, hx, eta,T,x0);





Wea1_XXssN_zL =IRFLow(:,1);
Wea2_XXssN_zL =IRFLow(:,2);
debtforhighincome_XXssN_zL  =IRFLow(:,3);
debtforlowincome_XXssN_zL  =IRFLow(:,4);
l1XXssN_zL  =IRFLow(:,5);
l2XXssN_zL  =IRFLow(:,6);
ProbDef1_ssN_zL  =IRFLow(:,7);
ProbDef2_ssN_zL =IRFLow(:,8);




sig = 0.01;

IRF_2ndLow = irf_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, T); 




Wea1_XXssN_z2L =IRF_2ndLow(:,1);
Wea2_XXssN_z2L =IRF_2ndLow(:,2);
debtforhighincome_XXssN_z2L  =IRF_2ndLow(:,3);
debtforlowincome_XXssN_z2L  =IRF_2ndLow(:,4);
l1XXssN_z2L  =IRF_2ndLow(:,5);
l2XXssN_z2L  =IRF_2ndLow(:,6);
ProbDef1_ssN_z2L  =IRF_2ndLow(:,7);
ProbDef2_ssN_z2L =IRF_2ndLow(:,8);




%% Model Simulation (low kappa)

sig = 0.01;
num_sim = 10000;
e = randn(num_sim, 1);

x0 = [0];
[YL, XL] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);

smplL = 0.9*num_sim:num_sim;


%% IRF





figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,1,1)
%plot(1:T,debtforhighincome_XXssN_z,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,debtforhighincome_XXssN_z2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,debtforhighincome_XXssN_zL,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,debtforhighincome_XXssN_z2L,'--' ,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Debt - High Income','FontSize',20)
legend({'High Collateral Req. kappa = 0.97','Low Collateral Req. kappa = 0.8'},'FontSize',20,'Location','best')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 

subplot(2,1,2)
%plot(1:T,debtforlowincome_XXssN_z,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,debtforlowincome_XXssN_z2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
%plot(1:T,debtforlowincome_XXssN_zL,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,debtforlowincome_XXssN_z2L,'--','LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Debt - Low Income','FontSize',20)
legend({'High Collateral Req. kappa = 0.97','Low Collateral Req. kappa = 0.8'},'FontSize',20,'Location','best')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 



%%






figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 

subplot(2,1,1)
plot(1:T,Wea1_XXssN_z2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,Wea1_XXssN_z2L,'--' ,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;

plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wealth - High Income','FontSize',20)
legend({'High Collateral Req. kappa = 0.97','Low Collateral Req. kappa = 0.8'},'FontSize',20,'Location','best')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 

subplot(2,1,2)
plot(1:T,Wea2_XXssN_z2,'LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,Wea2_XXssN_z2L,'--','LineWidth',3,'Color', [0, 158, 169]/255); hold on;
plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Wealth - Low Income','FontSize',20)
legend({'High Collateral Req. kappa = 0.97','Low Collateral Req. kappa = 0.8'},'FontSize',20,'Location','best')

axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 


%% Plot consumtion and saving response to investment schock
% 


% 0.09*(1+ratiorate_Rp2)*100
%  8.9880
%   11.7305
%   13.5029
%   14.5705
%   15.1324
%   15.3373
%   15.2956
%   15.0887
%   14.7762
%   14.4014
%   13.9952
%   13.5793
%   13.1688
%   12.7735
%   12.3997
%   12.0510
%   11.7292
%   11.4346
%   11.1667
%   10.9246
%   10.7065
%   10.5110
%   10.3362
%   10.1804
%   10.0417
%    9.9186
%    9.8095
%    9.7129
%    9.6275
%    9.5521
%    9.4855
%    9.4268
%    9.3751
%    9.3296
%    9.2895
%    9.2543
%    9.2233
%    9.1960
%    9.1721
%    9.1510



%% plot SIMULATE HIGH

 

%figure(10); 
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
subplot(2,1,1)

plot(smpl,Y(smpl,1),'LineWidth',1.5,'Color','#223F59'); hold on;
plot(smpl,Y(smpl,3),'LineWidth',1.5,'Color','#5890A6');
hold on;
legend({'Wealth ', 'Debt '},'FontSize',20,'Location','best','interpreter','latex')
%title('Output','FontSize',30);
ax = gca; % to enlarge the size of years and percentage on the x and y axis
ax.FontSize = 18;
ax.XGrid = 'off';
ax.YGrid = 'on';
ylabel('\% points','interpreter','latex');
xlabel('Time','FontSize',18,'interpreter','latex');
title('A. High Income Locations','FontSize',18,'interpreter','latex')

subplot(2,1,2)
plot(smpl,Y(smpl,2),'LineWidth',1.5,'Color','#223F59'); hold on;
plot(smpl,Y(smpl,4),'LineWidth',1.5,'Color','#5890A6') 
hold on;
legend('Wealth','Debt','FontSize',20,'Location','best','interpreter','latex')
ax = gca; % to enlarge the size of years and percentage on the x and y axis
ax.FontSize = 18;
ax.XGrid = 'off';
ax.YGrid = 'on';
ylabel('\% points','interpreter','latex');
xlabel('Time','FontSize',18,'interpreter','latex');
title('B. Low Income Locations','FontSize',18,'interpreter','latex')

export_fig PLotWealthDebt1 -m2.5 -painters


%%

figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

subplot(2,1,1)
plot(smpl,Y(smpl,1),'LineWidth',1.5,'Color','#223F59'); hold on;
plot(smplL,YL(smplL,1),'LineWidth',1.5,'Color','#5890A6');
hold on;
legend('High $\kappa$','Low $\kappa$','FontSize',20,'Location','best','interpreter','latex')
ax = gca; % to enlarge the size of years and percentage on the x and y axis
ax.FontSize = 18;
ax.XGrid = 'off';
ax.YGrid = 'on';
ylabel('\% points','interpreter','latex');
xlabel('Time','FontSize',18,'interpreter','latex');
title('A. Wealth in High Income Locations','FontSize',18,'interpreter','latex')

subplot(2,1,2)
plot(smpl,Y(smpl,2),'LineWidth',1.5,'Color','#223F59'); hold on ;
plot(smplL,YL(smplL,2),'LineWidth',1.5,'Color','#5890A6');
hold on;
legend('High $\kappa$','Low $\kappa$','FontSize',20,'Location','best','interpreter','latex')
ax = gca; % to enlarge the size of years and percentage on the x and y axis
ax.FontSize = 18;
ax.XGrid = 'off';
ax.YGrid = 'on';
ylabel('\% points','interpreter','latex');
xlabel('Time','FontSize',18,'interpreter','latex');
title('B. Wealth in Low Income Locations','FontSize',18,'interpreter','latex')


%saveas(gcf,'proxyPlot1','epsc2')

export_fig PLotWealth1 -m2.5 -painters






%%   Estimation in level (correlation between debt an delta wealth)

TotalWealth1meanLevel =416445.61;
TotalWealth2meanLevel =240529.20;

% this change in wealth(Delta W) not wealth 
wealth1mean = 106656.696776; %416445.61 (level); %PSID see code
wealth2mean =  56543.910020;%240529.20 (level) ;  %PSID see code

debt1mean = 15340.272153;%PSID data see code high income
debt2mean = 12370.732966 ;%PSID data see code high income

beta1_corre = 1.55; % see above 
beta2_corre = 0.81; % see above

%% mean of debt simulation when kappa is high   (from model)
Debt1_mean_HighKappa = mean(debt1mean *(1+Y(smpl,3)));
Debt2_mean_HighKappa = mean(debt2mean *(1+Y(smpl,4)));
%% mean of debt simulation when kappa is low
Debt1_mean_LowKappa  = mean(debt1mean *(1+YL(smplL,3)));
Debt2_mean_LowKappa  = mean(debt2mean *(1+YL(smplL,4)));

%% mean of  wealth 
TWealth1_mean_HighKappa = mean(TotalWealth1meanLevel *(1+Y(smpl,1)));
TWealth2_mean_HighKappa = mean(TotalWealth2meanLevel *(1+Y(smpl,2)));
%% mean of wealth simulation when kappa is low
TWealth1_mean_LowKappa  = mean(TotalWealth1meanLevel *(1+YL(smplL,1)));
TWealth2_mean_LowKappa  = mean(TotalWealth2meanLevel *(1+YL(smplL,2)));


%%


Diffweal1=(TWealth1_mean_HighKappa - TWealth1_mean_LowKappa);
Diffweal2=(TWealth2_mean_HighKappa - TWealth2_mean_LowKappa);

Diffweal1percent=(TWealth1_mean_HighKappa - TWealth1_mean_LowKappa)/TWealth1_mean_LowKappa*100;
Diffweal2percent=(TWealth2_mean_HighKappa - TWealth2_mean_LowKappa)/TWealth2_mean_LowKappa*100;


%%   Estimation in log (elasticity between debt an delta wealth)

LOGwealth1mean = 9.81; %PSID see ols regression
LOGwealth2mean =  10.42; %PSID see ols regression

beta1elas=0.25; % estimate of elasiticity see paper table results
beta2elas= 0.13; % estimate of elasticity
%% Predict \delta wealth (levels) using the simulation from the model 
 %% high \kappa
WWealth_Level_HighIncome_kappaHigh = wealth1mean + beta1_corre * Debt1_mean_HighKappa ; 
WWealth_Level_LowIncome_kappaHigh  = wealth2mean + beta2_corre * Debt2_mean_HighKappa;
 %% low \kappa
WWealth_Level_HighIncome_kappaLow = wealth1mean + beta1_corre * Debt1_mean_LowKappa; 
WWealth_Level_LowIncome_kappaLow  = wealth2mean + beta2_corre * Debt2_mean_LowKappa;
%%  
Diff1 = (WWealth_Level_HighIncome_kappaHigh - WWealth_Level_HighIncome_kappaLow) ;
Diff2 = (WWealth_Level_LowIncome_kappaHigh - WWealth_Level_LowIncome_kappaLow) ;

Diff1percent = (WWealth_Level_HighIncome_kappaHigh - WWealth_Level_HighIncome_kappaLow)/WWealth_Level_HighIncome_kappaLow*100 ;
Diff2percent = (WWealth_Level_LowIncome_kappaHigh - WWealth_Level_LowIncome_kappaLow)/WWealth_Level_LowIncome_kappaLow*100 ;


%% Predict Log delta wealth using the simulation from the model
 %% high \kappa
WWealth_Log_HighIncome_kappaHigh = LOGwealth1mean +  beta1elas *  log(Debt1_mean_HighKappa);
WWealth_Log_LowIncome_kappaHigh  = LOGwealth2mean +  beta2elas *  log(Debt2_mean_HighKappa);
 %% low \kappa
WWealth_Log_HighIncome_kappaLow  = LOGwealth1mean +  beta1elas *  log(Debt1_mean_LowKappa);
WWealth_Log_LowIncome_kappaLow   = LOGwealth2mean +  beta2elas *  log(Debt2_mean_LowKappa);



%% report table





if true 
fid = fopen('tabCounterfactual.xls','w+');
% shocks
fprintf(fid, '%30s \t %30s   \t %30s          \n','','High Income States','Low Income States');
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Level of Wealth (mean) ',TotalWealth1meanLevel  , TotalWealth2meanLevel );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Level of Change in Wealth (mean) ', wealth1mean , wealth2mean );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Log of Change in Wealth (mean)',  LOGwealth1mean , LOGwealth2mean);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Level of Debt (mean)', debt1mean , debt2mean);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Wealth Accumulation \Delta W ',0, 0);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Low Collateral Requirements \kappa', WWealth_Level_HighIncome_kappaLow , WWealth_Level_LowIncome_kappaLow );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','High Collateral Requirements \kappa',WWealth_Level_HighIncome_kappaHigh  , WWealth_Level_LowIncome_kappaHigh);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Difference ', Diff1 , Diff2 );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Difference in percent', Diff1percent , Diff2percent );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Wealth W (mean)',0, 0);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Low Collateral Requirements \kappa ',TWealth1_mean_LowKappa, TWealth2_mean_LowKappa);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','High Collateral Requirements \kappa ',TWealth1_mean_HighKappa, TWealth2_mean_HighKappa);
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Difference ', Diffweal1, Diffweal2 );
fprintf(fid, '%30s \t %12.4f \t %12.4f        \n','Difference in percent ', Diffweal1percent,Diffweal2percent );
fclose(fid);
end
 
 