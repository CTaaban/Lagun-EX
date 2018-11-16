%% Spectral wave analysis for OSSI (pressure sensor) time-series

close all; clear all; clc;

%% Summary:
%{
0. Supply [INPUT] parameters
1. Import
2. Pressure signal 
3. SSE
4. Wave Spectrum
5. Export wave characteristics 
%}

%% 0. Supply [INPUT] parameters
ossi = 'OSSI4'; % OSSI number 
type = '1';     % Type 1 for saturday-tuesday, 2 for wednesday-friday
log = '006';    % 006 for 25 sept
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

Patm=1.0378;    % [INPUT] Atmospheric pressure
rho=1027;       % [INPUT] density water 
zpt = 0.23;     % [INPUT] sensor elevation
g=9.81;         % [INPUT] gravitational acceleration
Fs=20;          % [INPUT] Sampling frequency

%% 1. Import:
if ossi == 'OSSI4' 
    if type == '1'
        data0 = csvread(['../Data/OSSI/',ossi,'/1st Measurement/','WLOG_',log,'.csv'],3);
    elseif type == '2'
        data0 = csvread(['../Data/OSSI/',ossi,'/2nd Measurement/','WLOG_',log,'.csv'],3);
    end
else
    data0 = csvread(['../Data/OSSI/',ossi,'/WLOG_',log,'.csv'],3);
end

data0 = data0(:,1:12); % For removing 13th column
data=data0';
b=numel(data0);
A=reshape(data,[b,1])+1;

addpath('../Functions/');

%% 2. Pressure signal:
figure();
plot(A)
title('Pressure sensor data OSSI4 (25/09/2018)')
xlabel('Sampling steps')
ylabel('Pressure [bar]')
set(gca,'fontsize',14)
%ylim([1.19,1.35])

if saveOn
    print(['../Data/Images/OSSI/',ossi,'_Pt.png'],'-dpng');
    close all;
end

%% 3. SSE:
dt=1/Fs; 
B=A-Patm;
pwave=detrend(B);
R_eta=pwave*10^5/(rho*g);

steps=round((20*60*60)*1/6);

t0 = '17:29'; % [INPUT] start time of time-series w.r.t. start of OSSI
tend = '17:36'; % [INPUT] end time of time-series w.r.t. start of OSSI

step1_HW = round((str2num(t0(1:2))+str2num(t0(4:5))/60)*length(A)/(length(A)/(Fs*3600)));
step2_HW = round((str2num(tend(1:2))+str2num(tend(4:5))/60)*length(A)/(length(A)/(Fs*3600)));

hHW=mean(B(step1_HW:step2_HW)*10^5/(rho*g));
eHW=p2sse(R_eta(step1_HW:step2_HW),dt,zpt,hHW);
etaHW=detrend(eHW);

t_step= 1/(20*60*60);
h_step1_HW = (step1_HW)/(20*60*60);
h_step2_HW = (step2_HW)/(20*60*60);
t_HW = [h_step1_HW:t_step:h_step2_HW];

figure(); hold on;
plot(t_HW,etaHW,'b','linewidth',1)
title(['Water surface elevation OSSI4 (25/09/2018: ',t0,'-',tend,')'])
xlabel('Hour of the day')
ylabel('\eta [m]') 
set(gca,'fontsize',14)
xlim([h_step1_HW,h_step2_HW])
a = (linspace(min(t_HW),max(t_HW),7));
xticks(a);
%ylim([0.35,2])

if saveOn
    print(['../Data/Images/OSSI/',ossi,'_Zt.png'],'-dpng');
    close all;
end

%% 4. Wave spectral analysis:
clc;
eta = etaHW; % [INPUT] time-series
p = 6; % [INPUT] number of blocks

n = length(eta); 
nfft = n/p;

%% 4.1 Number of blocks: 
%{
% 1 block:
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
% 4. Blocks vs. frequency resolution:
figure();     
subplot(2,1,1);
plot(info(:,1), info(:,3), 'r-', 'linewidth',1.5); hold on;
plot(info(:,1), info(:,4), 'g-', 'linewidth',1.5); hold on;
plot([0, p_i(end)], [1, 1], 'linewidth',1.5);
legend('ConfLow', 'ConfUpper');
ylabel('Confidence Intervals');
xlim([1, p_i(end)]);
title('Varying number of blocks');
set(gca,'fontsize',14);
subplot(2,1,2);
plot(info(:,1),info(:,2), 'linewidth',2);
xlabel('Number of blocks (p)'); 
ylabel('Frequency resolution (Hz)');
xlim([1, p_i(end)]);
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
hold on;
plot(f,E,'linewidth',2); 
xlabel('f [Hz]'); ylabel('E [m^2]'); 
xlim([0,1]);ylim([0,0.4]);
title(['Variance Density Spectrum (OSSI4)']); hold on;
plot(f,ConfUpper*E,'r--','linewidth',2);
plot(f,ConfLow*E,'g--','linewidth',2);
legend(['VDS (',num2str(p), ' blocks)'],['Upper limit 90% confidence interval'],['Lower limit 90% confidence interval']);
set(gca,'fontsize',14);

if saveOn
    print(['../Data/Images/OSSI/',ossi,'_Ef_1.png'],'-dpng');
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

%% 4.2 Check Gaussian distribution fit:
% define time vector [s]
time = (0:length(eta)-1)/Fs; 

% definition of the bins of elevation used to estimate the probability density function from the data
bin_edges = -1.5:0.1:1.5;              % edges of the bins
bin_width = bin_edges(2)-bin_edges(1); % width of the elevation bins

% we divide the surface elevation data in bins
bin_counts = histc(eta,bin_edges); % bin_counts is the number of elements in each of the pre-defined bins 
                                     
% we calculate the probability density function for the sea surface elevation

pdf_from_data = bin_counts/(sum(bin_counts)*bin_width); 

% This empirical probability density function will be compared with the gaussian distribution
m0n = var(eta);  % calculation zero-th order moment (=variance)
pdf_gaussian = 1/(2*pi*m0n)^0.5*exp(-bin_edges.^2/(2*m0n)); % Gaussian distribution

figure,subplot(2,1,1),plot(time,eta)
title(['Water surface elevation ', num2str(ossi), ' (25/09/2018: ',t0,'-',tend,')'])
xlabel('t [s]'); ylabel('\eta [m]')
xlim([0,time(end)]);ylim([-1,1]);
set(gca,'fontsize',14);
subplot(2,1,2), bar(bin_edges,pdf_from_data,'histc')
hold on,plot(bin_edges,pdf_gaussian,'r','linewidth',2)
xlabel('\eta [m]'); ylabel('p(\eta)')
xlim([-1.5,1.5]);
title('Probability distribution')
legend('Histogram time-series','Gaussian distribution')

set(gca,'fontsize',14)

if saveOn
    print(['../Data/Images/OSSI/',ossi,'_Pn.png'],'-dpng');
    close all;
end

Pr_emp = sum(bin_width * pdf_from_data);
Pr_gauss = trapz(bin_edges, pdf_gaussian);

disp(' ');
disp(['Position: ', num2str(ossi)]);
disp(['P_emp = ', num2str(Pr_emp)]);
disp(['P_gauss = ', num2str(Pr_gauss)]);

%% 4.3 Wave periods: 
F_nyq = Fs/2;
F_nyq = 1.526251526251526;
m2 = spectral_moment(f,E,0,F_nyq,2);
m1 = spectral_moment(f,E,0,F_nyq,1);
m0 = spectral_moment(f,E,0,F_nyq,0);
m_1 = spectral_moment(f,E,0,F_nyq,-1);
disp(' ');

cutoff = 1;
index_max = find(E==max(E(cutoff:end)));
f_p = f(index_max);
T_p = 1 / f_p;


T_m02 = sqrt(m0/m2);
T_m01 = m0/m1;
T_m_10 = m_1/m0;

disp(['T_m02 = ', num2str(T_m02), ' s']);
disp(['T_m01 = ', num2str(T_m01), ' s']);
disp(['T_m_10 = ', num2str(T_m_10), ' s']);
disp(['T_p = ', num2str(T_p), ' s']);

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

xlim([0 1.8]);ylim([0,0.4]);
xlabel('f [Hz]'); ylabel('E [m^2/Hz]')
legend('E','f_{m01}','f_{m02}', 'f_{m-10}')
title(['Variance Density Spectrum (OSSI4)']); hold on;
legend(['VDS (', num2str(p), ' blocks)'],['Upper limit 90% confidence interval'],['Lower limit 90% confidence interval'],'f_{m-10}', 'f_{m01}','f_{m02}');
set(gca,'fontsize',14);

if saveOn
    print(['../Data/Images/OSSI/',ossi,'_Ef_2.png'],'-dpng');
    close all;
end

%% 4.4 Wave heights:
m0  = spectral_moment(f,E,0,F_nyq,0);
Hm0 = 4 * sqrt(m0);
H_mean = sqrt(2*pi*m0);
H_rms = sqrt(8*m0);

disp(' ');
disp(['H_mean = ',num2str(H_mean),'m']);
disp(['H_rms = ',num2str(H_rms),'m']);
disp(['Hm0 = ',num2str(Hm0),'m']);
disp(' ');

%% 5. Export wave characteristics:
if saveOn
   E_2 = E;
   f_2 = f;
   ConfLow_2 = ConfLow;
   ConfUpper_2 = ConfUpper;
   T_p_2 = T_p;
   T_m02_2 = T_m02;
   T_m01_2 = T_m01;
   T_m_10_2 = T_m_10;
   m0_2 = m0;
   Hm0_2 = Hm0;
   H_mean_2 = H_mean;
   H_rms_2 = H_rms;
   save(['../Data/mat/Spectrum/',GPS_dataset(4:8),'/',map_dataset],'E_2','f_2',...
   'ConfLow_2','ConfUpper_2','T_p_2','T_m02_2','T_m01_2','T_m_10_2',...
   'm0_2','Hm0_2','H_mean_2','H_rms_2','-append');
end
