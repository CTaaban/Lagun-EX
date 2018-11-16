%% Transform spherical map (RTPA) to sliced map with grouped points for each LiDAR ray

% Author: Chahid Taaban
% Created: 23-10-2018

clc, clear all, close all;

count = 0; % used for not shifting GPS data multiple times
count2 = 0; % used for not subtracting LiDAR origin multiple times

%% Others:
%{
% Determining delta T:
delta = 1.0;
for i=2:length(RTP(:,2))
    if abs(RTP(i,2)-RTP(i-1,2)) < delta && abs(RTP(i,2)-RTP(i-1,2)) > 10^-5
        delta = abs(RTP(i,2)-RTP(i-1,2));
    end
end
disp(['Delta: ', num2str(delta)])
%}

%% Summary:
% 0. Supply [INPUT] parameters 
% 1. Import spherical dataset (RTPA)
% 2. Reduce wave data 
% 3. Make wave slices
% 4  Manipulate 3D matrix 
% 5. Export RTPA_XYZ_rays 

%% 0. Supply [INPUT] parameters
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

%% 1. Import spherical dataset (RTPA)
struct = load(['../Data/mat/RTPA/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
GPS = dlmread(['../Data/GPS/',GPS_dataset(4:8),'/',GPS_dataset,'_Adapted.txt'],'',1,0); 
scarp = dlmread(['../Data/GPS/Other/Scarp_xyz.txt'],',',1,0); 
zero_lvl = dlmread(['../Data/GPS/Other/Zero_lvl_xyz.txt'],',',1,0); 
RTPA = struct.RTPA;
clear struct;
addpath('../Functions/');

% Shift GPS coordinates according to cloudCompare values:
if count == 0
    shift = [72700.00,452600.00]; % [x-coordinate,y-coordinate]
    GPS(:,1:2) = GPS(:,1:2) - shift;
    scarp(:,1:2) = scarp(:,1:2) - shift;
    zero_lvl(:,1:2) = zero_lvl(:,1:2) - shift;
    count = -999;
end

ref_pole_A = GPS(1,1:3);
ref_pole_B = GPS(2,1:3);
ref_pole_C = GPS(3,1:3);
LiDAR = GPS(4,1:3);

% Make LiDAR position Origin
if count2 == 0
    scarp(:,1:2) = scarp(:,1:2) - LiDAR(1:2); 
    zero_lvl(:,1:2) = zero_lvl(:,1:2) - LiDAR(1:2); 
    ref_pole_A(:,1:2) = ref_pole_A(:,1:2) - LiDAR(1:2); 
    ref_pole_B(:,1:2) = ref_pole_B(:,1:2) - LiDAR(1:2); 
    ref_pole_C(:,1:2) = ref_pole_C(:,1:2) - LiDAR(1:2); 
    count2 = -999;
end

 %% 2. Reduce wave data
R_min_lim = 20;  % [INPUT] used for reducing data-set
R_max_lim = 60;  % [INPUT] used for reducing data-set
T_min_lim = -180; % [INPUT] used for reducing data-set
T_max_lim = +180; % [INPUT] used for reducing data-set

T_min = min(RTPA(:,2));
T_max = max(RTPA(:,2));

count3 = 1;
rng('default');
clear RTPA2;

for i=1:size(RTPA,1)
    if rand < 1.0 % to reduce randomly points
        if RTPA(i,1) >= R_min_lim && RTPA(i,1) <= R_max_lim % reduce based on R
            if RTPA(i,2) >= T_min_lim && RTPA(i,2) <= T_max_lim % reduce based on T
                RTPA2(count3,1) = RTPA(i,1);
                RTPA2(count3,2) = RTPA(i,2);
                RTPA2(count3,3) = RTPA(i,3);
                RTPA2(count3,4) = RTPA(i,4);
                count3 = count3 + 1;
            end
        end
    end
end

% Set RGB values based on reflectivity A:
RGB = zeros(length(RTPA2(:,4)),3);
%RGB(:,1) = RTPA2(:,4).^0.4;
RGB(:,2) = RTPA2(:,4).^0.15; 
RGB(:,3) = 1-RTPA2(:,4).^0.15; 

% Plot Reduced data-set:
figure()
polarscatter(0, 0, 150, 'r', 'filled'); hold on;
polarscatter(RTPA2(:,2)/180*pi, RTPA2(:,1), 5, RGB ,'filled'); hold on;
polarscatter(0, 0, 150, 'r', 'filled'); hold on;
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [round(T_min), round(T_max)];
ax.ThetaTick = round(linspace(round(T_min), round(T_max), 4));
ax.RLim = [0 R_max_lim+10];
ax.RTick = round(linspace(0, 100, 21));
R_num = 9;
ax.RTick = linspace(R_min_lim, R_max_lim, R_num);

set(gca, 'fontsize',14);
title(['2D polar view: ', map_dataset]);

R_labels = {};
for i=0:R_num
    R_labels{i+1} = [num2str(R_min_lim + i*(R_max_lim - R_min_lim)/(R_num-1)), 'm'];
end
rticklabels(R_labels)
ax.RColor = [0, 0.5, 0];

off = round((T_max-T_min)/4);
theta_labels = {};
for i=0:(T_max - T_min)/(off)+1
    theta_labels{i+1} = [num2str(round(T_min + i*(T_max - T_min)/(off))),char(176)];
end
thetaticklabels(theta_labels)
ax.ThetaColor = [0, 0.4, 0.45];
%ax.GridColor = 'red';

legend('LiDAR','Waves');
grid on;
if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_reduced.png'],'-dpng');
    close all;
end

%% 3.1 Make wave slices
T_range = T_max - T_min; % range of theta, not exaclty as MetaData value 
T_delta = MetaData{index,7}; % From MetaData from manual of LiDAR
%T_delta = 3; % used only for five slices img
rays = ceil((T_range/T_delta)+1); % Number of rays of LiDAR scan

% 3D matrix: all points (rows) with RTPA values (columns) of each wave slice (3rd dimension):
clear RTPA_rays;
RTPA_rays= zeros(size(RTPA2,1),4,rays); 

max_count_i = 1;

for i=1:size(RTPA2,1)
    if RTPA2(i,2)==T_max
        count_j = rays;
    else
        count_j = ceil((RTPA2(i,2)-T_max)/(T_min-T_max)*rays); % register to which ray a certain point belongs 
    end

    count_i = find(~RTPA_rays(:,1,count_j),1); % find empty row to fill with new point 
    RTPA_rays(count_i,1,count_j) = RTPA2(i,1); % Fill R-value
    RTPA_rays(count_i,2,count_j) = RTPA2(i,2); % Fill T-value
    RTPA_rays(count_i,3,count_j) = RTPA2(i,3); % Fill P-value
    RTPA_rays(count_i,4,count_j) = RTPA2(i,4); % Fill A-value
    
    if count_i > max_count_i 
        max_count_i = count_i; % register max. row number of total rays
    end
end

RTPA_rays(max_count_i+1:end,:,:) = []; % reducing empty rows of all rays equally, longest rows of ray is normative

%% 3.2 Plot Wave Rays1
mod_rays = 10; % [INPUT] used for reducing amount of rays plotted  

% Set RGB values based on reflectivity A:
RGB = zeros(max_count_i,3,size(RTPA_rays,3));
%RGB(:,1,:) = RTPA_rays(:,4,:).^0.4;
RGB(:,2,:) = RTPA_rays(:,4,:).^0.15; 
RGB(:,3,:) = 1-RTPA_rays(:,4,:).^0.15; 

figure()
polarscatter(0, 0, 150, 'r', 'filled'); hold on;

for i = 1:size(RTPA_rays,3)
    if mod(i,mod_rays)==0
        polarscatter(RTPA_rays(:,2,i)/180*pi, RTPA_rays(:,1,i), 5, RGB(:,:,i), 'filled'); hold on;
    end
end

polarscatter(0, 0, 150, 'r', 'filled'); hold on;
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [round(T_min), round(T_max)];
ax.ThetaTick = round(linspace(round(T_min), round(T_max), 4));
ax.RLim = [0 R_max_lim+10];
R_num = 9;
ax.RTick = linspace(R_min_lim, R_max_lim, R_num);

set(gca, 'fontsize',14);
title(['2D polar view: ', map_dataset]);

R_labels = {};
for i=0:R_num
    R_labels{i+1} = [num2str(R_min_lim + i*(R_max_lim - R_min_lim)/(R_num-1)), 'm'];
end
rticklabels(R_labels)
ax.RColor = [0, 0.5, 0];

off = round((T_max-T_min)/4);
theta_labels = {};
for i=0:(T_max - T_min)/(off)+1
    theta_labels{i+1} = [num2str(round(T_min + i*(T_max - T_min)/(off))), char(176)];
end
thetaticklabels(theta_labels)
ax.ThetaColor = [0, 0.4, 0.45];
%ax.GridColor = 'red';

legend('LiDAR','Waves');
grid on;

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_rays1.png'],'-dpng');
    close all;
end

%% Plot Wave Rays2 (colors)
colors = {'m','c','r','g','b','k','y'}; % [INPUT]
mod_rays = 1; % [INPUT] used for reducing amount of rays plotted  

figure()
polarscatter(0, 0, 150, 'r', 'filled'); hold on;
for i = 1:size(RTPA_rays,3)
    if mod(i,mod_rays)==0
            color_index = mod(i,size(colors,2));
        if color_index == 0
            color_index = size(colors,2); 
        end
        polarscatter(RTPA_rays(:,2,i)/180*pi, RTPA_rays(:,1,i), 5, colors{color_index} ,'filled'); hold on;
    end
end
polarscatter(0, 0, 150, 'r', 'filled');
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [round(T_min), round(T_max)];
ax.ThetaTick = round(linspace(round(T_min), round(T_max), 4));
ax.RLim = [0 R_max_lim+10];
R_num = 9;
ax.RTick = linspace(R_min_lim, R_max_lim, R_num);

set(gca, 'fontsize',14);
title(['2D polar view: ', map_dataset]);

R_labels = {};
for i=0:R_num
    R_labels{i+1} = [num2str(R_min_lim + i*(R_max_lim - R_min_lim)/(R_num-1)), 'm'];
end
rticklabels(R_labels)
ax.RColor = [0, 0.5, 0];

off = round((T_max-T_min)/4);
theta_labels = {};
for i=0:(T_max - T_min)/(off)+1
    theta_labels{i+1} = [num2str(round(T_min + i*(T_max - T_min)/(off))),char(176)];
end
thetaticklabels(theta_labels)
ax.ThetaColor = [0, 0.4, 0.45];
%ax.GridColor = 'red';

legend('LiDAR','Waves');
grid on;

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_rays2.png'],'-dpng');
    close all;
end

%% 4  Manipulate 3D matrix   

% Ascending rows:
Asc_var = 1; % [INPUT] - Options: 1=R, 2=T, 3=P, 4=A
for i=1:size(RTPA_rays,3)
    RTPA_rays_V2(:,:,i) = sortrows(RTPA_rays(:,:,i),Asc_var);
end

% From Spherical back to cartesian coordinates:
for i=1:size(RTPA_rays,3)
    %[x,y,z] = sph2cart(RTPA(:,2), RTPA(:,3),RTPA(:,1));

    x = RTPA_rays_V2(:,1,i).*cos(RTPA_rays_V2(:,2,i)*pi/180+pi).*sin(RTPA_rays_V2(:,3,i)*pi/180); 
    y = RTPA_rays_V2(:,1,i).*sin(RTPA_rays_V2(:,2,i)*pi/180+pi).*sin(RTPA_rays_V2(:,3,i)*pi/180);
    z = RTPA_rays_V2(:,1,i).*cos(RTPA_rays_V2(:,3,i)*pi/180);

    XYZ(:,1:3,i) = [x,y,z];
end

% Combine both 3D matrices to one 
RTPAXYZ_rays(:,1:4,:) = RTPA_rays_V2;
RTPAXYZ_rays(:,5:7,:) = XYZ;

% Set zero values to NaN
RTPAXYZ_rays(RTPAXYZ_rays == 0) = NaN;

%% 4.2 2D plot of Rays3 (XYZ):
mod_rays = 10; % [INPUT] to reduce number of rays plotted

figure();
plot(0,0,'r.','markersize',50); hold on;
plot(scarp(:,1), scarp(:,2),'-','color',[0.6350, 0.0780, 0.1840],'linewidth',1.5); hold on;
plot(zero_lvl(:,1),zero_lvl(:,2),'-','color',[0, 0, 1],'linewidth',1.5); hold on;
plot(ref_pole_A(1),ref_pole_A(2),'.','color',[0.3010, 0.7450, 0.9330],'markersize',40); hold on;
plot(ref_pole_B(1),ref_pole_B(2),'.','color',[0, 0.4470, 0.7410],'markersize',40); hold on;
plot(ref_pole_C(1),ref_pole_C(2),'.','color',[0.4660, 0.6740, 0.1880],'markersize',40); hold on;

for i = 1:size(XYZ,3)
    if mod(i,mod_rays)==0
       color_index = mod(i,size(colors,2));
    if color_index == 0
       color_index = size(colors,2); 
    end
    scatter(XYZ(:,1,i),XYZ(:,2,i),5, colors{color_index}); hold on;
    end
end

plot(0,0,'r.','markersize',50); hold on;

set(gca, 'fontsize',14);
title(['2D Topdown view: ', map_dataset, ' (Rays)']);
%xlim([min(X)-5,0]);ylim([0,max(Y)+5]);
xlabel('X [m]');ylabel('Y [m]');
legend('LiDAR','Scarp','0m NAP','Ref. A','Ref. B','Ref. C','Waves');
grid on;

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_rays3.png'],'-dpng');
    close all;
end

%% 5. Export RTPA_XYZ
if saveOn
    save(['../Data/mat/RTPAXYZ_rays/',GPS_dataset(4:8),'/',map_dataset],'RTPAXYZ_rays');
end

