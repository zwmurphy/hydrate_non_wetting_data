clear all; clc; close all

t=xlsread('HWRP05_intial_conditions.xlsx','Initial_saturation','E2:E674');
water=xlsread('HWRP05_intial_conditions.xlsx','Initial_saturation','F2:F674');
gas=xlsread('HWRP05_intial_conditions.xlsx','Initial_saturation','G2:G674');

x=[158 158];
y=[0 60];
%% Plot vs time
hold on

plot(t,water,'-b','Linewidth',1.25)
plot(t,gas,'-r','Linewidth',1.25)
plot(x,y,'--k','Linewidth',1.25)
xlabel('time (hrs)'); ylabel('volume (ml)')
legend('Water removed','Gas injected','location','best')
xlim([0 inf]); %ylim([0 1]);

set(gca, ...
  'FontName','Arial',...
  'FontSize',20,...
  'Box'         , 'off'     , ...
  ...%'xticklabel'  , {[]}    ,...
  'TickDir'     , 'out'     , ...
  ...%'TickLength'  , [.02,.02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XGrid'       , 'off'     ,...
  'YGrid'       , 'off'      , ...
  'XColor'      , 'k', ...
  'YColor'      , 'k', ...
  ...%'YTick'       , 0:500:2500, ...
  'LineWidth'   , 1);