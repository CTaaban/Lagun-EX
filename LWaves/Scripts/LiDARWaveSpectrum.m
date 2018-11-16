%% Spectral wave analysis for LiDAR time-series (RZt)
%
% Author: Chahid Taaban
% Created: 25-10-2018

clear all; close all; clc;

%% Summary:
%{
0. Supply [INPUT] parameters 
1. Import dataset (RZt) 
2. Selecting time-series
3. Interpolate missing values
4. De-trend time-series 
5. Wave envelope 
6. Wave-by-wave analysis 
7. Wave spectral analysis
8. Export wave characteristics
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
struct = load(['../Data/mat/RZt/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
R = struct.R_struc;
Z = struct.Z_struc;
t = struct.t_struc;
t_0 = struct.t_0;
clear struct;
addpath('../Functions/');

%% 2. Selecting time-series:
R_set = 45; % [INPUT] radius of time-series
R_index = find(R==R_set);
Z_1 = Z(R_index,:);

%% 3. Interpolate missing values:
a = nnz(isnan(Z_1));
Z_2 = fillmissing(Z_1,'linear');
b = nnz(isnan(Z_2));

%% 4. De-trend time-series:
Z_3 = detrend(Z_2);
idx_first = find(sum(~isnan(Z_1),1) > 0, 1 ,'first');
idx_last = find(sum(~isnan(Z_1),1) > 0, 1 , 'last');

Z_4 = zeros(1,length(Z_3));
t_4 = zeros(1,length(Z_3));
t_0_V2 = t_0 + idx_first*t_delta;

for i=1:length(Z_3)
    if i >= idx_first && i <= idx_last
        t_4(i) = t(i);
        Z_4(i) = Z_3(i); 
    else
        t_4(i) = nan;
        Z_4(i) = nan;
    end
end

t_5 = t(Z_4 < nanmedian(Z_4) + 3*nanstd(Z_4));
Z_5 = Z_4(Z_4 < nanmedian(Z_4) + 3*nanstd(Z_4));

figure('units','normalized','outerposition',[0 0 1 1]);
%figure()
plot(t,Z_1,'k','linewidth',1.5); hold on;
plot(t,Z_2,'r','linewidth',1.5); hold on;
plot(t,Z_3,'b','linewidth',1.5); hold on;
plot(t,Z_4,'g','linewidth',1.5); hold on;
plot(t_5,Z_5,'m','linewidth',1.5); hold on;
plot([0, t(end)], [0, 0],'r', 'linewidth',2);

plot(t,Z_1,'k','linewidth',1.5); hold on;
%ylim([z_min,z_max]);
set(gca,'fontsize',24);
grid;
title(['Wave time-series at R = ', num2str(R_set),'m (', map_dataset,')']);
ylabel('\eta [m]');
xlabel('Time [s]');
legend(['Original'],['Interpolated'],['De-trended'],['Selected '],['Outliers removed']);

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_Zt_V3.png'],'-dpng');
    close all;
end

t_4(isnan(t_4)) = []; 
Z_4(isnan(Z_4)) = []; 

%% 5. Wave envelope: 
Z_5 = Z_4(Z_4 > nanmedian(Z_4) + 1*nanstd(Z_4));
eta = Z_5;

H = max(eta)-min(eta);
disp(['Wave envelope: ', num2str(H) ,'m']);

%% 6. Wave-by-wave analysis:
eta = Z_4;
[H_ind, T_ind] = zero_crossing(eta,1/t_delta);
T_0 = mean(T_ind);
H_mean = mean(H_ind);
H_s = significant_wave_height(H_ind);
H_rms = rms_wave_height(H_ind);
H_max = max(H_ind);

disp(' ');
disp(['T_0 = ', num2str(T_0),'s']);
disp(['H_mean = ', num2str(H_mean),'m']);
disp(['H_rms = ', num2str(H_rms),'m']);
disp(['H_s = ', num2str(H_s),'m']);
disp(['H_max = ', num2str(H_max), ' m']);
disp(' ');

%% 7. Wave spectral analysis:
close all; clc;
eta = Z_3; % [INPUT] time-series
p = 6; % [INPUT] number of blocks

n = length(eta); 
nfft = n/p;
Fs = 1/t_delta;

%% 7.1 Number of blocks: 
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
[E,f,ConfLow,ConfUpper] = wave_spectrum(eta,nfft,Fs);
f_delta = f(2,1) - f(1,1);
t_delta = 1/Fs;
f_nyquist = 1 / (2 * t_delta);

disp(['f_delta = ', num2str(f_delta), ' Hz']);
disp(['f_nyquist = ', num2str(f_nyquist), ' Hz']);

figure(); 
plot(f,E,'linewidth',2); 
xlabel('f [Hz]'); ylabel('E [m^2]'); 
xlim([0,1]);ylim([0,0.4]);
title(['Variance Density Spectrum (LiDAR)']); hold on;
plot(f,ConfUpper*E,'r--','linewidth',2);
plot(f,ConfLow*E,'g--','linewidth',2);
legend(['VDS (', num2str(p), ' blocks)'],['Upper limit 90% confidence interval'],['Lower limit 90% confidence interval']);
set(gca,'fontsize',14);

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_Ef_1.png'],'-dpng');
    close all;
end

%{ 
% Total energy:
rho = 1000;
variance = trapz(f, E);
Energy = variance*9.81*rho;
disp(' ');
disp(['E = ', num2str(Energy), ' J/m^2 or kg/s^2']);
%}

%% 7.2 Check Gaussian distribution fit:
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

%% 7.3 Wave periods: 
F_nyq = Fs/2;
m2 = spectral_moment(f,E,0,F_nyq,2);
m1 = spectral_moment(f,E,0,F_nyq,1);
m0 = spectral_moment(f,E,0,F_nyq,0);
m_1 = spectral_moment(f,E,0,F_nyq,-1);
disp(' ');

E_max = max(E);
i = 1;
for j = 1:length(E)
    if E(j)==E_max
        break;
    end
    i = i + 1;
end
T_p1 = 1 / f(i);

index_max = find(E==max(E(5 :end)));
f_p = f(index_max);
T_p2 = 1 / f_p;
T_m02 = sqrt(m0/m2);
T_m01 = m0/m1;
T_m_10 = m_1/m0;

disp(['T_m02 = ', num2str(T_m02), ' s']);
disp(['T_m01 = ', num2str(T_m01), ' s']);
disp(['T_m_10 = ', num2str(T_m_10), ' s']);
disp(['T_p = ', num2str(T_p1), ' s']);
%disp(['T_p2 = ', num2str(T_p2), ' s']);

f01 = T_m01^-1;
f02 = T_m02^-1;
fminus10 = T_m_10^-1;

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

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_Ef_2.png'],'-dpng');
    close all;
end

%% 7.4 Wave heights: 
m0  = spectral_moment(f,E,0,F_nyq,0);
Hm0 = 4 * sqrt(m0);
H_mean = sqrt(2*pi*m0);
H_rms = sqrt(8*m0);

disp(' ');
disp(['H_mean = ',num2str(H_mean),'m']);
disp(['H_rms = ',num2str(H_rms),'m']);
disp(['Hm0 = ',num2str(Hm0),'m']);
disp(' ');

%% 8. Export time-series:
if saveOn
   E_1 = E;
   f_1 = f;
   ConfLow_1 = ConfLow;
   ConfUpper_1 = ConfUpper;
   T_p1_1 = T_p1;
   T_p2_1 = T_p2;
   T_m02_1 = T_m02;
   T_m01_1 = T_m01;
   T_m_10_1 = T_m_10;
   m0_1 = m0;
   Hm0_1 = Hm0;
   H_mean_1 = H_mean;
   H_rms_1 = H_rms;
   save(['../Data/mat/Spectrum/',GPS_dataset(4:8),'/',map_dataset],'E_1','f_1',...
   'ConfLow_1','ConfUpper_1','T_p1_1','T_p2_1','T_m02_1','T_m01_1','T_m_10_1',...
   'm0_1','Hm0_1','H_mean_1','H_rms_1','-append');
save(['../Data/mat/Precision/',GPS_dataset(4:8),'/',map_dataset],'eta');
end


