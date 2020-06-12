clear
%% during hydrate formation
t=xlsread('Mass Balance.xlsx','A2:A2853');
T=xlsread('Mass Balance.xlsx','C2:C2853');
P=xlsread('Mass Balance.xlsx','F2:F2853');
P_up=xlsread('Mass Balance.xlsx','E2:E2853');
P_down=xlsread('Mass Balance.xlsx','F2:F2853');
Vpump=xlsread('Mass Balance.xlsx','D2:D2853');
save('HWRP05.mat','t','T','P','P_up','P_down','Vpump')
%% during gas removal
% t2=xlsread('Mass Balance.xlsx','A835:A896');
% T2=xlsread('Mass Balance.xlsx','C835:C896');
% P2=xlsread('Mass Balance.xlsx','B835:B896');
% delta_w2=xlsread('Mass Balance.xlsx','J835:J896');
% save('HWRP01.mat','t2','T2','P2','delta_w2')