clear
%% during hydrate formation
t=xlsread('Mass Balance.xlsx','A2:A2658');
T=xlsread('Mass Balance.xlsx','C2:C2658');
P=xlsread('Mass Balance.xlsx','B2:B2658');
Vpump=xlsread('Mass Balance.xlsx','D2:D2658');
save('HWRP07.mat','t','T','P','Vpump')

% t2=xlsread('Mass Balance.xlsx','A835:A896');
% T2=xlsread('Mass Balance.xlsx','C835:C896');
% P2=xlsread('Mass Balance.xlsx','B835:B896');
% delta_w2=xlsread('Mass Balance.xlsx','J835:J896');
% save('HWRP01.mat','t2','T2','P2','delta_w2')