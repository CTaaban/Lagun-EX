%% Sensitivity analysis of precision on spectral wave characteristics of LiDAR (Precision)

% Author: Chahid Taaban
% Created: 15-11-2018

clear all; close all; clc;

%{
TODO:
Calcualte Noise and NSR auto
%}

%% Summary:
%{
0. Supply [INPUT] parameters 
1. Import dataset  
2. Wave spectral analysis
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

% RUN ID, Date, Time, Resolution, angular range, angular increment, scan, duration, angular increment
MetaData = {'LW007', '25-09-2018', '17:18', 3.1, 10, '01:32', 0.018;
    'LW008', '25-09-2018', '17:29', 1.6, 10, '06:04', 0.009;
    'LW010', '26-09-2018', '16:46', 3.1, 10, '01:32', 0.018;
    'LW011', '26-09-2018', '16:51', 1.6,  4, '02:27', 0.009;
    'LW017', '27-09-2018', '17:23', 3.1, 10, '01:32', 0.018;
    'LW018', '27-09-2018', '17:27', 1.6, 10, '06:04', 0.009;
    'LW020', '27-09-2018', '17:37', 0.8,  4, '04:51', 0.004;
    'LW021', '27-09-2018', '18:01', 1.6,  4, '02:27', 0.009;
    'LW022', '27-09-2018', '18:21', 3.1, 40, '06:01', 0.018;
    'LW023', '27-09-2018', '18:29', 3.1, 10, '01:32', 0.018;
    }; % [INPUT]

% Find row index for dataset of input
for i=1:size(MetaData,1)
    if MetaData{i,1} == map_dataset
        index = i;
        break;
    else
        index = -999;
    end
end

t_dur = str2double(MetaData{index,6}(1:2))*60+str2double(MetaData{index,6}(4:5)); % Total scan duration in [s]
t_delta = t_dur*MetaData{index,7}/MetaData{index,5}; % time between two rays

%% 1. Import data: 
struct = load(['../Data/mat/Precision/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
struct2 = load(['../Data/mat/Spectrum/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
eta = struct.eta;
clear struct;
addpath('../Functions/');

% OSSI4 wave characteristics: 
E_2 = struct2.E_2;
f_2 = struct2.f_2;
ConfLow_2 = struct2.ConfLow_2;
ConfUpper_2 = struct2.ConfUpper_2;
T_p_2 = struct2.T_p_2;
T_m02_2 = struct2.T_m02_2;
T_m01_2 = struct2.T_m01_2;
T_m_10_2 = struct2.T_m_10_2;
m0_2 = struct2.m0_2;
Hm0_2 = struct2.Hm0_2;
H_mean_2 = struct2.H_mean_2;
H_rms_2 = struct2.H_rms_2;

%% 2. Wave spectral analysis:
deltaZ = 0.0036597686; % [INPUT] Precision error in [m]

rng('default');

a = -1;
b = 1;
r = ((b-a).*rand(1,length(eta)) + a)*deltaZ;

eta_1 = eta;
eta_2 = eta_1 + r;

p = 6; % [INPUT] number of blocks
n = length(eta); 
nfft = n/p;
Fs = 1/t_delta;

% 7.1 Number of blocks: 
%{
% 1 Block:
p_0 = 1;
nfft_0 = n/p_0;
[E,f,ConfLow,ConfUpper] = wave_spectrum(eta,nfft_0,Fs);
figure(); 
plot(f,E,'linewidth',1); hold on;
plot(f,ConfUpper*E,'r--','linewidth',1);hold on;
plot(f,ConfLow*E,'g--','linewidth',1); hold on;
xlabel('f (Hz)'); ylabel('E (m^2)'); 
legend(['Raw estimate'],['Upper confidence limit (5%)'],['Lower confidence limit (5%)']);
title('Variance Density Spectrum'); hold on;
set(gca,'fontsize',14);
%}

i = 0;
info = zeros(5, 5);

for p_i=[1,2,3,4,5,6,10,12]  
    i = i + 1;
    info(i,1) = p_i;
    
    nfft_i = n/p_i;
    [E,f,ConfLow,ConfUpper] = wave_spectrum(eta,nfft_i,Fs);
    
    f_delta = f(2,1) - f(1,1);
    info(i,2) = f_delta;
    
    %{
    % Plots:
    figure; 
    plot(f,E,'linewidth',1.5); 
    xlabel('f (Hz)'); ylabel('E (m^2)'); 
    title(['Estimate Variance Density Spectrum (', num2str(p), ' blocks)']); hold on;
    
    plot(f,ConfUpper*E,'r--','linewidth',1.5);
    plot(f,ConfLow*E,'g--','linewidth',1.5);
    set(gca,'fontsize',14);
    %}
    info(i,3) = ConfLow;
    info(i,4) = ConfUpper;
    
    variance = trapz(f, E);
    info(i,5) = variance;
    
end


%{
% Variance:
p=[1,2,3,4,5,6,10,12];  
disp(' ');
var_eta = std(eta)^2;
disp(['var(eta) = ', num2str(var_eta), ' m^2']);
for i=1:length(info)
    disp(['Variance (',num2str(p(i)),' blocks) = ',num2str(info(i,5)), ' m^2']);
end
%}

%{
% Blocks vs. frequency resolution:
figure();     
subplot(2,1,1);
plot(info(:,1), info(:,3), 'r-', 'linewidth',1.5); hold on;
plot(info(:,1), info(:,4), 'g-', 'linewidth',1.5); hold on;
plot([0, p(end)], [1, 1], 'linewidth',1.5);
legend('ConfLow', 'ConfUpper');
ylabel('Confidence Intervals');
xlim([1, p(end)]);
title('Varying number of blocks');
set(gca,'fontsize',14);
subplot(2,1,2);
plot(info(:,1),info(:,2), 'linewidth',2);
xlabel('Number of blocks (p)'); 
ylabel('Frequency resolution (Hz)');
xlim([1, p(end)]);
set(gca,'fontsize',14);
%}

% 5. Final block with Tp, delta_f and f_nyq:
[E_1,f_1,ConfLow_1,ConfUpper_1] = wave_spectrum(eta_1,nfft,Fs);
[E_2,f_2,ConfLow_2,ConfUpper_2] = wave_spectrum(eta_2,nfft,Fs);

f_delta = f_1(2,1) - f_1(1,1);
t_delta = 1/Fs;
f_nyquist = 1 / (2 * t_delta);

disp(['f_delta = ', num2str(f_delta), ' Hz']);
disp(['f_nyquist = ', num2str(f_nyquist), ' Hz']);

figure(); hold on;
plot(f_1,E_1,'linewidth',2); 
plot(f_2,E_2,'linewidth',2); 
xlabel('f [Hz]'); ylabel('E [m^2]'); 
%xlim([0,1]);ylim([0,0.4]);
title(['Sensititivy of variance Density Spectrum (LiDAR)']); hold on;
legend(['VDS (', num2str(p), ' blocks)'],['VDS with elevation increase/decrease']);
set(gca,'fontsize',14);

if saveOn
    print(['../Data/Images/Precision/Ef_4.png'],'-dpng');
    %close all;
end

%{ 
% Total energy:
rho = 1000;
variance = trapz(f, E);
Energy = variance*9.81*rho;
disp(' ');
disp(['E = ', num2str(Energy), ' J/m^2 or kg/s^2']);
%}

%% 2.2 Check Gaussian distribution fit:
%{
% define time vector [s]
time = (0:length(eta)-1)/Fs; 

% other useful variables
radius = 45; % radius of the time-series

% definition of the bins of elevation used to estimate the probability density function from the data
bin_edges = -1.5:0.1:1.5;              % edges of the bins
bin_width = bin_edges(2)-bin_edges(1); % width of the elevation bins

% we divide the surface elevation data in bins
bin_counts = histc(eta,bin_edges); % bin_counts is the number of elements in each of the pre-defined bins 
                                     
% we calculate the probability density function for the sea surface elevation
pdf_from_data = bin_counts/(sum(bin_counts)*bin_width); 

% This empirical probability density function will be compared with  the gaussian distribution
m0n = var(eta);  % calculation zero-th order moment (=variance)
pdf_gaussian = 1/(2*pi*m0n)^0.5*exp(-bin_edges.^2/(2*m0n)); % Gaussian distribution

t_end_min = str2num(MetaData{index,3}(4:5))+str2num(MetaData{index,6}(1:2));
figure,subplot(2,1,1),plot(time,eta)
title(['Water surface elevation at R = ' num2str(radius) 'm (',MetaData{index,2},': ',MetaData{index,3},'-',MetaData{index,3}(1:2),':',num2str(t_end_min),')'])
xlabel('t [s]'); ylabel('\eta [m]')
xlim([0,375]);ylim([-1,1]);
set(gca,'fontsize',14);
subplot(2,1,2), bar(bin_edges,pdf_from_data,'histc')
hold on,plot(bin_edges,pdf_gaussian,'r','linewidth',2)
xlabel('\eta [m]'); ylabel('p(\eta)')
xlim([-1.5,1.5]);
title('Probability distribution')
legend('Histogram time-series','Gaussian distribution')
set(gca,'fontsize',14)

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_Pn.png'],'-dpng');
    close all;
end

Pr_emp = sum(bin_width * pdf_from_data);
Pr_gauss = trapz(bin_edges, pdf_gaussian);

disp(' ');
disp(['Position: R=', num2str(radius),'m']);
disp(['P_emp = ', num2str(Pr_emp)]);
disp(['P_gauss = ', num2str(Pr_gauss)]);
%}

%% 2.3 Wave periods: 
F_nyq = Fs/2;

m2_1 = spectral_moment(f_1,E_1,0,F_nyq,2);
m1_1 = spectral_moment(f_1,E_1,0,F_nyq,1);
m0_1 = spectral_moment(f_1,E_1,0,F_nyq,0);
m_1_1 = spectral_moment(f_1,E_1,0,F_nyq,-1);

m2_2 = spectral_moment(f_2,E_2,0,F_nyq,2);
m1_2 = spectral_moment(f_2,E_2,0,F_nyq,1);
m0_2 = spectral_moment(f_2,E_2,0,F_nyq,0);
m_1_2 = spectral_moment(f_2,E_2,0,F_nyq,-1);

index_max_1 = find(E_1==max(E_1(1:end)));
f_p_1 = f_1(index_max_1);
T_p1_1 = 1/f_p_1;
T_m02_1 = sqrt(m0_1/m2_1);
T_m01_1 = m0_1/m1_1;
T_m_10_1 = m_1_1/m0_1;

index_max_2 = find(E_2==max(E_2(1:end)));
f_p_2 = f_2(index_max_2);
T_p1_2 = 1/f_p_2;
T_m02_2 = sqrt(m0_2/m2_2);
T_m01_2 = m0_2/m1_2;
T_m_10_2 = m_1_2/m0_2;

disp(' ');
disp('Delta z = 0:')
disp(['T_m02 = ', num2str(T_m02_1), ' s']);
disp(['T_m01 = ', num2str(T_m01_1), ' s']);
disp(['T_m_10 = ', num2str(T_m_10_1), ' s']);
disp(['T_p = ', num2str(T_p1_1), ' s']);
disp(' ');
disp(['Delta z = ', num2str(deltaZ),':'])
disp(['T_m02 = ', num2str(T_m02_2), ' s']);
disp(['T_m01 = ', num2str(T_m01_2), ' s']);
disp(['T_m_10 = ', num2str(T_m_10_2), ' s']);
disp(['T_p = ', num2str(T_p1_2), ' s']);
disp(' ');


delta_T_m02 = max([T_m02_1-T_m02_2]);
delta_T_m01 = max([T_m01_1-T_m01_2]);
delta_T_m_10 = max([T_m_10_1-T_m_10_2]);
delta_T_p = max([T_p1_1-T_p1_2]);

disp(' ');
disp('Delta T:');
disp(['T_m02 = ', num2str(delta_T_m02), ' s']);
disp(['T_m01 = ', num2str(delta_T_m01), ' s']);
disp(['T_m_10 = ', num2str(delta_T_m_10), ' s']);
disp(['T_p = ', num2str(delta_T_p), ' s']);

per_T_m02 = delta_T_m02/T_m02_1*100;
per_T_m01 = delta_T_m01/T_m01_1*100;
per_T_m_10 = delta_T_m_10/T_m_10_1*100;
per_T_p = delta_T_p/T_p1_1*100;

disp(' ');
disp('Percentage T:');
disp(['T_m02 = ', num2str(per_T_m02),'%']);
disp(['T_m01 = ', num2str(per_T_m01),'%']);
disp(['T_m_10 = ', num2str(per_T_m_10),'%']);
disp(['T_p = ', num2str(per_T_p),'%']);


f01_1 = T_m01_1^-1;
f02_1 = T_m02_1^-1;
fminus10_1 = T_m_10_1^-1;

f01_2 = T_m01_2^-1;
f02_2 = T_m02_2^-1;
fminus10_2 = T_m_10_2^-1;

%{
figure();
plot(f,E,'linewidth',2); hold on;
plot(f,ConfUpper*E,'r--', 'linewidth',2); hold on;
plot(f,ConfLow*E,'g--', 'linewidth',2); hold on;
plot([fminus10 fminus10],[0 1],'m', 'linewidth',2); hold on;
plot([f01 f01],[0 1],'k', 'linewidth',2); hold on;
plot([f02 f02],[0 1],'r', 'linewidth',2); hold on;

xlim([0 1]);ylim([0,0.4]);
xlabel('f [Hz]'); ylabel('E [m^2/Hz]')
legend('E','f_{m01}','f_{m02}', 'f_{m-10}')
title(['Variance Density Spectrum (LiDAR)']); hold on;
legend(['VDS (', num2str(p), ' blocks)'],['Upper limit 90% confidence interval'],['Lower limit 90% confidence interval'],'f_{m-10}', 'f_{m01}','f_{m02}');
set(gca,'fontsize',14);

%}

if saveOn
    print(['../Data/Images/Precision/Ef_5.png'],'-dpng');
    close all;
end

%% 2.4 Wave heights: 

m0_1  = spectral_moment(f_1,E_1,0,F_nyq,0);
Hm0_1 = 4 * sqrt(m0_1);
H_mean_1 = sqrt(2*pi*m0_1);
H_rms_1 = sqrt(8*m0_1);

m0_2  = spectral_moment(f_2,E_2,0,F_nyq,0);
Hm0_2 = 4 * sqrt(m0_2);
H_mean_2 = sqrt(2*pi*m0_2);
H_rms_2 = sqrt(8*m0_2);


disp(' ');
disp(['Delta z = 0:'])
disp(['H_mean = ',num2str(H_mean_1),'m']);
disp(['H_rms = ',num2str(H_rms_1),'m']);
disp(['Hm0 = ',num2str(Hm0_1),'m']);
disp(' ');
disp(['Delta z = +', num2str(deltaZ),':'])
disp(['H_mean = ',num2str(H_mean_2),'m']);
disp(['H_rms = ',num2str(H_rms_2),'m']);
disp(['Hm0 = ',num2str(Hm0_2),'m']);


delta_H_mean = max([H_mean_1-H_mean_2]);
delta_H_rms = max([H_rms_1-H_rms_2]);
delta_H_m0 = max([Hm0_1-Hm0_2]);

disp(' ')
disp('Delta H:')
disp(['H_mean = ',num2str(delta_H_mean),'m']);
disp(['H_rms = ',num2str(delta_H_rms),'m']);
disp(['Hm0 = ',num2str(delta_H_m0),'m']);

per_H_mean = delta_H_mean/H_mean_1*100;
per_H_rms = delta_H_rms/H_rms_1*100;
per_H_m0 = delta_H_m0/Hm0_1*100;

disp(' ')
disp('Percentage H:')
disp(['H_mean = ',num2str(per_H_mean),'%']);
disp(['H_rms = ',num2str(per_H_rms),'%']);
disp(['Hm0 = ',num2str(per_H_m0),'%']);

