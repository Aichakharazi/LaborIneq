load('macrotimeseries1.mat')
%%
load('Recessions.mat')
start=Recessions(:,1);
finish=Recessions(:,2);
%%

startdate = datenum('Q1-1964','QQ-yyyy');%startdate = datenum('04-1997','mm-yyyy');
enddate = datenum('Q4-2021','QQ-yyyy');%startdate = datenum('04-1997','mm-yyyy'); 
dt = linspace(startdate,enddate,232);

%%
figure('color','w') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 
[h1]=plot(dt,labobs);
shadeyears(start,finish);
hold on;
% Create a plot with 2 y axes using the plotyy function
[ax, h1, h2] = plotyy(dt,labobs,dt, yobs, 'plot');
axis tight;  
datetick('x','yyyy','keepticks');
% Add title and x axis label
xlabel('Years','FontSize',19);
set(ax(1),'XLim',[datenum(1964,3,1) datenum(2021,12,1)],'FontSize',15);  
set(ax(2),'Xlim',[datenum(1964,3,1) datenum(2021,12,1)],'FontSize',15); 
set(ax(1),'ycolor',[165,3,3]/255)
set(ax(2),'ycolor',[1,39,89]/255)
grid on
% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'),'String','Billions of Hours','FontSize',19,'Color',[165,3,3]/255 );
set(get(ax(2), 'Ylabel'),'String','Billions of Dollars','FontSize',19,'Color',[1,39,89]/255); %U.S. Dollars per Barrell
yt=get(gca,'YTick');
set(ax(1),'YTickLabel',sprintf('%.0f\n',yt))
set(ax(2),'YTickLabel',sprintf('%.0f\n',yt))
% set color of lines
set(h1,'LineWidth',2.5,'Color',[165,3,3]/255)
set(h2,'LineWidth',2.5,'Color',[1,39,89]/255)
lgd = legend([h1;h2],'Total Hours Worked','Gross Domestic Product','Location','best');
lgd.FontSize = 16;

%%
% Hours of Wage and Salary Workers on Nonfarm Payrolls: Total (TOTLQ)
% Billions of Hours,
% Seasonally Adjusted Annual Rate

% gdp Billions of Dollars, Seasonally Adjusted Annual Rate
