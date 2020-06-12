clc; close all
%% mass Balance for a hydrate system - Zach Murphy, Summer 2019
% 3 phases: l (liquid), g (vapor), h (solid)
% 3 components: w (water), m (methane), s (salt)

%% Assumtptions
% 1) Solubility of methane in water is very low (~.2%) and is ignored
% 2) Assume sample reaches T2 instantly
% 3) No water/salt in gas phase
% 4) No salt in hydrate phase
%% Import Data: MUST UPDATE FOR EACH EXPERIMENT!
load('HWRP07.mat');

%% Core parameters: Inputs %% Confirm for every experiment!!!!!!!
L=0.5588; %sample length (m)
r=0.01905; %sample radius (m)
A=pi*r^2; %sample area (m^2)
V_sample=A*L; %sample volume (m^3)
porosity=0.215; %sample porosity (-)
Vp=V_sample*porosity; %sample pore volume
%T=[23,6.5]; %temp (C) T1=start T2=temp for rest of experiment (assume to be instant/constant)
V_line=2; %volume of line btw pump and core (ml)
R=8.314; %gas constant (Pa m^3)/(mol k)
leak_rate=.9; %(ml/hr)
%% Known constraints: Constants (Can update/change if desired)
Mm=.01604; %MM of CH4 (kg/mol)
Xwg=0; %mass fraction of water in gas phase (assume to be zero)
Xml=0; %mass fraction of gas in liquid phase (assume to be zero)
Xsg=0; Xsh=0; %mass fraction of salt in gas and hydrate phases (assume to be zero)
Xmh=0.134; %mass fraction of gas in hydrate
Xwh=1-Xmh; %mass fraction of water in hydrate
Xmg=1; %mass fraction of methane in gas(gas is pure methane)
rho_h=912; %density of hydrate (kg/m^3)

%% Correct Nm for PT and pump compression
deltaNs=zeros(length(Vpump),1);
deltaNw=zeros(length(Vpump),1);
deltaNm=zeros(length(Vpump),1);

for ii=1:length(Vpump)
Vpump(ii)=(Vpump(ii)+leak_rate*t(ii));
deltaNm(ii)=((((Vpump(1)+V_line)*10^-6*P(1)*6894.75)/((T(1)+273.15)*R))-(((Vpump(ii)+V_line)*10^-6*P(ii)*6894.75)/((T(2)+273.15)*R)))*Mm;  % mass in gas (kg)
end


%% Experimental initial conditions (MUST KNOW THESE)!!!!!!!!!!
Sli=0.81; %inital brine saturation (know from mass balance/initial saturation procedure)
Sgi=0.19; %initial gas saturation (know from mass balance/initial saturation procedure)
Shi=1-Sli-Sgi; %initial hydrate saturation (Should be 0)
Xsli=0.04; Xwli=1-Xsli; % initial salinity of brine

%% initial conditions
rho_l_i=brine_density(P(1),T(1),Xsli/0.05844/(1-Xsli),Xml);  % initial brine density (kg/m^3)
rho_gi=CH4_EOS(T(1),P(1))*Mm; %initial gas density (kg/m^3)
Nwi=(Xwli*Sli*rho_l_i+Xwg*Sgi*rho_gi+Xwh*Shi*rho_h)*Vp; %mass of water in system (kg)
Nmi=(Xml*Sli*rho_l_i+Xmg*Sgi*rho_gi+Xmh*Shi*rho_h)*Vp; %mass of methane in system (kg)
Nsi=(Xsli*Sli*rho_l_i+Xsg*Sgi*rho_gi+Xsh*Shi*rho_h)*Vp; %mass of salt in system (kg)
%% Preallocate!
Ns=zeros(length(Vpump),1);
Nw=zeros(length(Vpump),1);
Nm=zeros(length(Vpump),1);
rho_g=zeros(length(Vpump),1);
rho_l=zeros(length(Vpump),1);
Sg=zeros(length(Vpump),1);
Sh=zeros(length(Vpump),1);
Sl=zeros(length(Vpump),1);
Xsl=zeros(length(Vpump),1);

%% Mass at each time step
for i=1:length(deltaNm)
 Nm(i)=Nmi+deltaNm(i); %mass CH4 (kg)
 Ns(i)=Nsi+deltaNs(i); %mass salt (kg)
 Nw(i)=Nwi+deltaNw(i); % mass water (kg)
end

%% Mass balance during hydrate formation (from solver)
for z=1:length(Nm)
    rho_g(z)=CH4_EOS(T(z),P(z))*Mm; %gas density at each step (kg/m^3)
    if z==1
        [Sg(z), Sh(z), Sl(z), Xsl(z)] = myfun(Nm(z),Ns(z),Nw(z),rho_l_i,rho_h,rho_g(z),Xml,Xwg,Xmg,Xsg,Xwh,Xmh,Xsh,Vp);
    else
        rho_l(z)=brine_density(P(z),T(z),Xsl(z-1)/0.05844/(1-Xsl(z-1)),Xml);  % new brine density (kg/m^3)
        [Sg(z), Sh(z), Sl(z), Xsl(z)] = myfun(Nm(z),Ns(z),Nw(z),rho_l(z),rho_h,rho_g(z),Xml,Xwg,Xmg,Xsg,Xwh,Xmh,Xsh,Vp);
    end
end

% %% Mass balance during gas removal -- Working on this....?? 
% % How can we determine mass balance if we do not know Nm, Ns, Nw? Can we force certain variables?
%  
% rho_g2=rho_g(end);
% for iz=1:length(delta_w2)
%     Nw2(iz)=Nwi+delta_w2(iz)/1000;
%     Sg2(iz)=0; %% Sl and Sh spike if this is immediately 0, can we figure out a way to measure? 
%     Xsl_2(iz)=.035; Xwl_2(iz)=1-Xsl_2(iz); %% all gas is removed to hydrate  / salinity is forced to .035
% 
% Nm2(iz)=(Nw2(iz)*Xmh*rho_h - Nw2(iz)*Xml*rho_l_i - Vp*Xmh*Xwl_2(iz)*rho_h*rho_l_i + Vp*Xml*Xwh*rho_h*rho_l_i + Sg2(iz)*Vp*Xmg*Xwh*rho_g2*rho_h - ...
%     Sg2(iz)*Vp*Xmh*Xwg*rho_g2*rho_h - Sg2(iz)*Vp*Xmg*Xwl_2(iz)*rho_g2*rho_l_i + Sg2(iz)*Vp*Xml*Xwg*rho_g2*rho_l_i + Sg2(iz)*Vp*Xmh*Xwl_2(iz)...
%     *rho_h*rho_l_i - Sg2(iz)*Vp*Xml*Xwh*rho_h*rho_l_i)/(Xwh*rho_h - Xwl_2(iz)*rho_l_i);
%  
%  
% Sh2(iz) =(Nw2(iz) - Vp*Xwl_2(iz)*rho_l_i - Sg2(iz)*Vp*Xwg*rho_g2 + Sg2(iz)*Vp*Xwl_2(iz)*rho_l_i)/(Vp*Xwh*rho_h - Vp*Xwl_2(iz)*rho_l_i);
%  
%  
% Sl2(iz) =-(Nw2(iz) - Vp*Xwh*rho_h - Sg2(iz)*Vp*Xwg*rho_g2 + Sg2(iz)*Vp*Xwh*rho_h)/(Vp*Xwh*rho_h - Vp*Xwl_2(iz)*rho_l_i);
% end
% 
%% Plot vs time
hold on

yyaxis left
plot(t,Sg,'-r','Linewidth',1.25)
plot(t,Sl,'-b','Linewidth',1.25)
plot(t,Sh,'-g','Linewidth',1.25)
% plot(t2,Sg2,'-r','Linewidth',1.25)
% plot(t2,Sl2,'-b','Linewidth',1.25)
% plot(t2,Sh2,'-g','Linewidth',1.25)
xlabel('time (hrs)'); ylabel('Saturation')
xlim([0.5 inf]); ylim([0 1]);

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

yyaxis right
plot(t,Xsl,'--k','Linewidth',1.25)
%plot(t2,Xsl_2,'--k','Linewidth',1.25)
ylabel('Salnity')
ylim([0,1]);
legend('Sg','Sl','Sh','Xsl','Location','best')
hold off

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