%% Converts Cartesian to Spherical Coordinates of 3D map
%
% Author: Chahid Taaban
% Created: 16-10-2018 

clc, clear all, close all;
    
count = 0; % used for not shifting GPS data multiple times 
count2 = 0; % used for not subtracting LiDAR origin multiple times

%% Summary:
% 0. Supply [INPUT] paramters
% 1. Import data: map & GPS
% 2. Subtract LiDAR position (x,y) from map
% 3. Plot cartesian dataset
% 4. Transform cartesian map coordinates to spherical 
% 5. Plot spherical dataset
% 6. Export spherical dataset 

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

%% 1. Import Data:
GPS = dlmread(['../Data/GPS/',GPS_dataset(4:8),'/',GPS_dataset,'_Adapted.txt'],'',1,0); 
map = dlmread(['../Data/txt/',GPS_dataset(4:8),'/',map_dataset, '.txt'], '',2,0);
scarp = dlmread(['../Data/GPS/Other/Scarp_xyz.txt'],',',1,0); 
zero_lvl = dlmread(['../Data/GPS/Other/Zero_lvl_xyz.txt'],',',1,0); 

% Shift GPS coordinates according to cloudCompare values:
if count == 0
    shift = [72700.00,452600.00]; % [x-coordinate,y-coordinate]
    GPS(:,1:2) = GPS(:,1:2) - shift;
    scarp(:,1:2) = scarp(:,1:2) - shift;
    zero_lvl(:,1:2) = zero_lvl(:,1:2) - shift;
    count = -999;
end

addpath('../Functions/');

%% 2. Subtract LiDAR position (only x,y) from map
ref_pole_A = GPS(1,1:3);
ref_pole_B = GPS(2,1:3);
ref_pole_C = GPS(3,1:3);
LiDAR = GPS(4,1:3);

if count2 == 0
    map(:,1:2) = map(:,1:2) - LiDAR(1:2);
    scarp(:,1:2) = scarp(:,1:2) - LiDAR(1:2); 
    zero_lvl(:,1:2) = zero_lvl(:,1:2) - LiDAR(1:2); 
    ref_pole_A(:,1:2) = ref_pole_A(:,1:2) - LiDAR(1:2); 
    ref_pole_B(:,1:2) = ref_pole_B(:,1:2) - LiDAR(1:2); 
    ref_pole_C(:,1:2) = ref_pole_C(:,1:2) - LiDAR(1:2); 
    count2 = -999;
end

%% 3. Plot cartesian dataset
X = map(:,1); % X-coordinate 
Y = map(:,2); % Y-coordinate 
Z = map(:,3); % Z-coordinate 
A = map(:,4); % Reflectivity scalar: Large R -> High reflectivity (green)  
A = (A - min(A))/(max(A)-min(A)); % Normalise A

% Define a color gradient based on reflectivity parameter used in plots
RGB = zeros(length(A),3);
%RGB(:,1) = A.^0.4;
RGB(:,2) = A.^0.15; 
RGB(:,3) = 1-A.^0.15;  

%% 3.1 2D plot:
figure();
scatter(X,Y,2,RGB,'filled'); hold on;
plot(scarp(:,1), scarp(:,2),'-','color',[0.6350, 0.0780, 0.1840],'linewidth',1.5); hold on;
plot(zero_lvl(:,1),zero_lvl(:,2),'-','color',[0, 0, 1],'linewidth',1.5); hold on;
plot(0,0,'r.','markersize',50); hold on;
plot(ref_pole_A(1),ref_pole_A(2),'.','color',[0.3010, 0.7450, 0.9330],'markersize',40); hold on;
plot(ref_pole_B(1),ref_pole_B(2),'.','color',[0, 0.4470, 0.7410],'markersize',40); hold on;
plot(ref_pole_C(1),ref_pole_C(2),'.','color',[0.4660, 0.6740, 0.1880],'markersize',40); hold on;
set(gca, 'fontsize',14);
title(['2D Topdown view: ', map_dataset]);
%xlim([min(X)-5,0]);ylim([0,max(Y)+5]);
xlabel('X [m]');ylabel('Y [m]');
legend('Waves','Scarp','0m NAP','LiDAR','Ref. A','Ref. B','Ref. C');
grid on;

%% 3.2 3D plot:
figure();
scatter3(X,Y,Z,2,RGB,'filled'); hold on;
plot(0,0,'r.','markersize',50);
set(gca, 'fontsize',14);
title(['3D view: ', map_dataset]);
xlim([min(X)-5,0]);ylim([0,max(Y)+5]);
xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]');
legend('Waves','LiDAR');
grid on;

%% 4. Transform cartesian map coordinates to spherical 
% Let r = radius (R), theta = angle on x-y plane (T), & phi = angle off z -axis (P). Then u get:

R = sqrt(X.^2 + Y.^2 + Z.^2);
T = atan(Y./X)/pi*180;
P = acos(Z./R)/pi*180; %degrees , gradient; 
RTPA = [R,T,P,A];

%% 5. Plot spherical dataset
if saveOn
    close all;
end

% Theta and Radius range for plot:
offset = 2; % [INPUT]
roundTo = 5; % [INPUT]
roundoff = round((max(T)-min(T))/roundTo);
theta_min = round((min(T)-offset)/roundoff)*roundoff; % [INPUT] use 300 for wave and 280 for runup view
theta_max = round((max(T)+offset)/roundoff)*roundoff; % [INPUT] use 320 for wave and 320 for runup view
off = round((theta_max-theta_min)/4);
R_min = 0; % [INPUT] 
R_max = max(round(R/5)*5); % [INPUT] 

%{
count3 = 1;
clear T2 R2 RGB2;
for i=1:length(T)
    if rand < 0.1
        if R(i) > 30 && R(i) < 31 && T(i) < -48 && T(i) > -49
        T2(count3,1) = T(i);
        R2(count3,1) = R(i);
        RGB2(count3,1) = RGB(i);
        count3 = count3 + 1;
        end
    end
end
%}

%% 5.1 2D view:
figure()
polarscatter(T/180*pi, R, 5, RGB, 'filled'); hold on;
%polarscatter(T2/180*pi, R2, 5, RGB2, 'filled'); hold on;
polarscatter(0, 0, 150, 'r', 'filled');
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaLim = [theta_min, theta_max];
ax.ThetaTick = round(linspace(theta_min, theta_max, (theta_max - theta_min)/(off)+1));
ax.RLim = [R_min R_max];
R_num = 16;
ax.RTick = linspace(R_min, R_max, R_num);

set(gca, 'fontsize',14);
title(['2D polar view: ', map_dataset]);

R_labels = {};
for i=0:R_num
    R_labels{i+1} = [num2str(R_min + i*(R_max - R_min)/(R_num-1)), 'm'];
end
rticklabels(R_labels)
ax.RColor = [0, 0.5, 0];

theta_labels = {};
for i=0:(theta_max - theta_min)/(off)+1
    theta_labels{i+1} = [num2str(theta_min + i*(theta_max - theta_min)/(off)),char(176)];
end
thetaticklabels(theta_labels);
ax.ThetaColor = [0, 0.4, 0.45];
%ax.GridColor = 'red';

legend('Waves','LiDAR');
grid on;
if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'.png'],'-dpng');
    close all;
end

%% 5.2 3D view:
%figure();
%polarplot3d(TRP);

%% 6. Export spherical dataset
if saveOn
    save(['../Data/mat/RTPA/',GPS_dataset(4:8),'/',map_dataset],'RTPA');
end

