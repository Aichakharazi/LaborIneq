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
 

y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_with_covidshock = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     % deviation from steady state  


//set_param_value('ddeltaCov',0); //
//stoch_simul(order=1,irf=11,nograph,nomoments,nocorr,nofunctions);
        
// set productivity shocks only e1 ;
shock_matrix(:,strmatch('e1',M_.exo_names,'exact')) = innovations.e1';
shock_matrix(:,strmatch('e2',M_.exo_names,'exact')) = zeros(1,n_points);
 

y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
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
plot(1:T,LM_hfV_exp1,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,LM_hfV_exp2,'LineWidth',3,'Color', '#0F5C8C'); hold on;
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
plot(1:T,LM_hnV_exp1,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,LM_hnV_exp2,'LineWidth',3,'Color', '#0F5C8C'); hold on;
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
plot(1:T,LM_lfV_exp1,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,LM_lfV_exp2,'LineWidth',3,'Color', '#0F5C8C'); hold on;
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
plot(1:T,LM_lnV_exp1,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,LM_lnV_exp2,'LineWidth',3,'Color', '#0F5C8C'); hold on;
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
plot(1:T,y_V_exp1,"-.",'LineWidth',3,'Color', '#0F95D7'); hold on;
plot(1:T,y_V_exp2,'LineWidth',3,'Color', '#0F5C8C'); hold on;
%plot(1:T,zeros(size(1:T)),'--' ,'LineWidth',3,'Color',[255, 0, 1]/255 ); 
%ylabel('% deviation from ss','FontSize',20)
xlabel('Time','FontSize',20);
title('Output','FontSize',20)
axis tight;  
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = 16; 
legend({'With Covid Shock','Without Covid Shock'},'FontSize',20,'Location','best')
