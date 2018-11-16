%% Validation of spectral wave characteristics of LiDAR w.r.t. OSSI's (Spectrum)

close all; clear all; clc;

%% Summary:
%{
0. Supply [INPUT] parameters
1. Import
2. Wave spectrum 
3. Wave frequencies
4. Wave periods 
5. Wave heights
6. Combined wave characteristics 
%}

%% 0. Supply [INPUT] parameters:
map_dataset = 'LW008'; % [INPUT] = name of 3D map dataset, LW020 densest so far and LW018 best
switch map_dataset % [INPUT] = name of GPS datasets (e.g. DAY_5_LW)
    case {'LW007', 'LW008'}
        GPS_dataset = 'LW_DAY_3'; 
    case {'LW010', 'LW011'}
        GPS_dataset = 'LW_DAY_4';
    case {'LW017', 'LW018', 'LW020', 'LW021', 'LW022', 'LW023'}
        GPS_dataset = 'LW_DAY_5';
end
saveOn = false; % [INPUT] boolean for saving img' and data

%% 1. Import:
addpath('../Functions/');
struct = load(['../Data/mat/Spectrum/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
% LiDAR:
E_1 = struct.E_1;
f_1 = struct.f_1;
ConfLow_1 = struct.ConfLow_1;
ConfUpper_1 = struct.ConfUpper_1;
T_p1_1 = struct.T_p1_1;
T_p2_1 = struct.T_p2_1;
T_m02_1 = struct.T_m02_1;
T_m01_1 = struct.T_m01_1;
T_m_10_1 = struct.T_m_10_1;
m0_1 = struct.m0_1;
Hm0_1 = struct.Hm0_1;
H_mean_1 = struct.H_mean_1;
H_rms_1 = struct.H_rms_1;

% OSSI4:
E_2 = struct.E_2;
f_2 = struct.f_2;
ConfLow_2 = struct.ConfLow_2;
ConfUpper_2 = struct.ConfUpper_2;
T_p_2 = struct.T_p_2;
T_m02_2 = struct.T_m02_2;
T_m01_2 = struct.T_m01_2;
T_m_10_2 = struct.T_m_10_2;
m0_2 = struct.m0_2;
Hm0_2 = struct.Hm0_2;
H_mean_2 = struct.H_mean_2;
H_rms_2 = struct.H_rms_2;
clear struct;

%% 2.1 Wave spectrum:

p = 6; % blocks
w1 = 5; % width plots

figure('units','normalized','outerposition',[0 0 1 1]);
%figure()
hold on;
plot(f_1,E_1,'r-','linewidth',w1); 
plot(f_2,E_2,'b-','linewidth',w1); 
xlabel('f [Hz]'); ylabel('E [m^2]'); 
xlim([0,1]);ylim([0,0.2]);
title(['Variance Density Spectra'])
legend(['LiDAR VDS (',num2str(p), ' blocks)'],['OSSI4 VDS (',num2str(p), ' blocks)'],'NumColumns',1);
set(gca,'fontsize',22);

if saveOn
    print(['../Data/Images/Validation/','Ef_0.png'],'-dpng');
    close all;
end

%% 2.2 Wave spectrum with confidence intervals:
p = 5; % blocks
w1 = 3; % width plots
w2 = 3; % width plots

figure('units','normalized','outerposition',[0 0 1 1]);
%figure()
hold on;
plot(f_1,E_1,'b-','linewidth',w1); 
plot(f_1,ConfUpper_1*E_1,'r-','linewidth',w1);
plot(f_1,ConfLow_1*E_1,'g-','linewidth',w1);

plot(f_2,E_2,'b-','linewidth',w2); 
plot(f_2,ConfUpper_2*E_2,'r-.','linewidth',w2);
plot(f_2,ConfLow_2*E_2,'g-.','linewidth',w2);
xlabel('f [Hz]'); ylabel('E [m^2]'); 
xlim([0,1]);%ylim([0,0.3]);
title(['Variance Density Spectra'])
legend(['LiDAR VDS (',num2str(p), ' blocks)'],'LiDAR upper limit 90% confidence interval','LiDAR lower limit 90% confidence interval',...
       ['OSSI4 VDS (',num2str(p), ' blocks)'],'OSSI4 upper limit 90% confidence interval','OSSI4 lower limit 90% confidence interval',...
       'NumColumns',2);
set(gca,'fontsize',22);

if saveOn
    print(['../Data/Images/Validation/','Ef_1.png'],'-dpng');
    close all;
end

%% 3. Wave frequencies:
b=[0, 0.4470, 0.7410];	          
o=[0.8500, 0.3250, 0.0980];	          
y=[0.9290, 0.6940, 0.1250];	          
pu=[0.4940, 0.1840, 0.5560];	          	
g=[0.4660, 0.6740, 0.1880];	      
c=[0.3010, 0.7450, 0.9330];	          
m=[0.6350, 0.0780, 0.1840];
w1 = 4;
w2 = 6;

f01_1 = T_m01_1^-1;
f02_1 = T_m02_1^-1;
fminus10_1 = T_m_10_1^-1;
fp1_1 = 1/T_p1_1;
fp2_1 = 1/T_p2_1;

f01_2 = T_m01_2^-1;
f02_2 = T_m02_2^-1;
fminus10_2 = T_m_10_2^-1;
fp_2 = 1/T_p_2;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(f_1,E_1,'r-','linewidth',w2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',0.01); 

plot([fp1_1 fp1_1],[0 1],'b-.','color',pu, 'linewidth',w1); 
plot([fminus10_1 fminus10_1],[0 1],'m-.','color',y, 'linewidth',w1);
plot([f01_1 f01_1],[0 1],'k-.','color',o, 'linewidth',w1); 
plot([f02_1 f02_1],[0 1],'r-.','color',b, 'linewidth',w1); 

plot(f_2,E_2,'b-','linewidth',w2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',0.01)
plot([fp_2 fp_2],[0 1],'g-','color',pu, 'linewidth',w1); 
plot([fminus10_2 fminus10_2],[0 1],'m-','color',y, 'linewidth',w1);
plot([f01_2 f01_2],[0 1],'k-','color',o, 'linewidth',w1); 
plot([f02_2 f02_2],[0 1],'r-','color',b, 'linewidth',w1); 

plot(f_1,E_1,'r-','linewidth',w2); 
plot(f_2,E_2,'b-','linewidth',w2); 

xlim([0 1]);ylim([0,0.2]);
xlabel('f [Hz]'); ylabel('E [m^2/Hz]')
legend('E','f_{m01}','f_{m02}', 'f_{m-10}')
title(['Comparison Variance Density Spectra']); hold on;
legend(['LiDAR VDS (',num2str(p), ' blocks)'],'LiDAR f_{p}','LiDAR f_{m-10}', 'LiDAR f_{m01}','LiDAR f_{m02}',...
       ['OSSI4 VDS (',num2str(p), ' blocks)'],'OSSI4 f_{p}','OSSI4 f_{m-10}', 'OSSI4 f_{m01}','OSSI4 f_{m02}',...
       'NumColumns',2);
set(gca,'fontsize',22);

if saveOn
    print(['../Data/Images/Validation/','Ef_2.png'],'-dpng');
    close all;
end

%% 4. Wave periods:
figure(); hold on
D = 300;

scatter(T_m02_1,T_m02_2,D,'filled');
scatter(T_m01_1,T_m01_2,D, 'filled');
scatter(T_m_10_1,T_m_10_2,D,'filled');
scatter(T_p1_1,T_p_2,D,'filled');
%scatter(T_p2_1,T_p_2,D,'filled');
plot([0,40], [0,40],'r-','linewidth',2)
xlabel('LiDAR wave periods [s]');
ylabel('OSSI4 wave periods [s]');
%xlim([0,20]);ylim([0,20])
title('Comparison of characteristic wave periods')
legend('T_{m02}','T_{m01}','T_{m-10}','T_{p}')
set(gca,'fontsize',14)
grid;

if saveOn
    print(['../Data/Images/Validation/','T.png'],'-dpng');
    close all;
end


disp(' ')
disp('LiDAR:')
disp(['T_m02_1 = ', num2str(T_m02_1), ' s']);
disp(['T_m01_1 = ', num2str(T_m01_1), ' s']);
disp(['T_m_10_1 = ', num2str(T_m_10_1), ' s']);
disp(['T_p_1 = ', num2str(T_p1_1), ' s']);
%disp(['T_p2_1 = ', num2str(T_p2_1), ' s']);
disp('OSSI4:')
disp(['T_m02_2 = ', num2str(T_m02_2), ' s']);
disp(['T_m01_2 = ', num2str(T_m01_2), ' s']);
disp(['T_m_10_2 = ', num2str(T_m_10_2), ' s']);
disp(['T_p_2 = ', num2str(T_p_2), ' s']);
disp('Delta T:')
disp(['T_m02 = ', num2str(T_m02_1-T_m02_2), ' s']);
disp(['T_m01 = ', num2str(T_m01_1-T_m01_2), ' s']);
disp(['T_m_10 = ', num2str(T_m_10_1-T_m_10_2), ' s']);
disp(['T_p = ', num2str(T_p1_1-T_p_2), ' s']);
disp('Percentage T:')
disp(['T_m02 = ', num2str((T_m02_1-T_m02_2)/T_m02_1*100), '%']);
disp(['T_m01 = ', num2str((T_m01_1-T_m01_2)/T_m01_1*100), '%']);
disp(['T_m_10 = ', num2str((T_m_10_1-T_m_10_2)/T_m_10_1*100), '%']);
disp(['T_p = ', num2str((T_p1_1-T_p_2)/T_p1_1*100), '%']);

%% 5. Wave heights:
figure(); hold on
D = 300;

scatter(H_mean_1,H_mean_2,D,'filled');
scatter(H_rms_1,H_rms_2,D,'filled');
scatter(Hm0_1,Hm0_2,D,'filled');
plot([0,1], [0,1],'r-','linewidth',2)
xlabel('LiDAR wave heights [m]');
ylabel('OSSI4 wave heights [m]');
title('Comparison of characteristic wave heights')
legend('H_{mean}','H_{rms}','H_{m0}')
set(gca,'fontsize',14)
grid;

if saveOn
    print(['../Data/Images/Validation/','H.png'],'-dpng');
    close all;
end

disp(' ');
disp('LiDAR:')
disp(['H_mean_1 = ',num2str(H_mean_1),'m']);
disp(['H_rms_1 = ',num2str(H_rms_1),'m']);
disp(['Hm0_1 = ',num2str(Hm0_1),'m']);
disp('OSSI4:')
disp(['H_mean_2 = ',num2str(H_mean_2),'m']);
disp(['H_rms_2 = ',num2str(H_rms_2),'m']);
disp(['Hm0_2 = ',num2str(Hm0_2),'m']);
disp('Delta H')
disp(['H_mean = ',num2str(H_mean_1-H_mean_2),'m']);
disp(['H_rms = ',num2str(H_rms_1-H_rms_2),'m']);
disp(['Hm0 = ',num2str(Hm0_1-Hm0_2),'m']);
disp('Percentage H')
disp(['H_mean = ',num2str((H_mean_1-H_mean_2)/H_mean_1*100),'%']);
disp(['H_rms = ',num2str((H_rms_1-H_rms_2)/H_rms_1*100),'%']);
disp(['Hm0 = ',num2str((Hm0_1-Hm0_2)/Hm0_1*100),'%']);

%% 6. Combined wave characteristics:
figure(); hold on
D = 300;

T_max = max([T_m02_1, T_m01_1, T_m_10_1, T_p1_1, T_m02_2, T_m01_2,T_m_10_2, T_p_2]);
H_max = max([H_mean_1, H_rms_1, Hm0_1, H_mean_2, H_rms_2, Hm0_2]);


scatter(T_m02_1/T_max,T_m02_2/T_max,D,'filled');
scatter(T_m01_1/T_max,T_m01_2/T_max,D, 'filled');
scatter(T_m_10_1/T_max,T_m_10_2/T_max,D,'filled');
scatter(T_p1_1/T_max,T_p_2/T_max,D,'filled');
scatter(H_mean_1/H_max,H_mean_2/H_max,D,'filled');
scatter(H_rms_1/H_max,H_rms_2/H_max,D,'filled');
scatter(Hm0_1/H_max,Hm0_2/H_max,D,'filled');
%scatter(T_p2_1,T_p_2,D,'filled');
plot([0,1.2], [0,1.2],'r-','linewidth',2)

errorbar(T_m02_1/T_max,T_m02_2/T_max,(T_m02_1-T_m02_2)/T_max,'horizontal','k-','color', b ,'linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
errorbar(T_m01_1/T_max,T_m01_2/T_max,(T_m01_1-T_m01_2)/T_max,'horizontal','k-','color', o ,'linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',10);
errorbar(T_m_10_1/T_max,T_m_10_2/T_max,(T_m_10_1-T_m_10_2)/T_max,'horizontal','k-','color', y ,'linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',15);
errorbar(T_p1_1/T_max,T_p_2/T_max,(T_p1_1-T_p_2)/T_max,'horizontal','k-','linewidth',2,'MarkerEdgeColor','k','color', pu,'MarkerFaceColor','k','CapSize',15);

errorbar(H_mean_1/H_max,H_mean_2/H_max,(H_mean_1-H_mean_2)/H_max,'horizontal','k-','linewidth',2,'MarkerEdgeColor','k','color', g,'MarkerFaceColor','k','CapSize',15);
errorbar(H_rms_1/H_max,H_rms_2/H_max,(H_rms_1-H_rms_2)/H_max,'horizontal','k-','linewidth',2,'MarkerEdgeColor','k','color', c,'MarkerFaceColor','k','CapSize',15);
errorbar(Hm0_1/H_max,Hm0_2/H_max,(Hm0_1-Hm0_2)/H_max,'horizontal','k-','linewidth',2,'MarkerEdgeColor','k','color', m,'MarkerFaceColor','k','CapSize',15);


xlabel('Normalised LiDAR values [-]');
ylabel('Normalsied OSSI4 values [-]');
%xlim([0,20]);ylim([0,20])
title('Comparison of wave characteristics')
legend('T_{m02}','T_{m01}','T_{m-10}','T_{p}','H_{mean}','H_{rms}','H_{m0}','NumColumns',2,'location','best')
set(gca,'fontsize',14)
grid;

if saveOn
    print(['../Data/Images/Validation/','HT.png'],'-dpng');
    close all;
end
