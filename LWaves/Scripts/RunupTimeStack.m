  %% Transform LiDAR 3D map to Wave Run-up Time-stack
%
% Author: Chahid Taaban
% Created: 25-10-2018

clear all; close all; clc;
count  = 0; % used for apply noise reduction only once

%% Notes for Khaled:

% Wave runup / waterline algorthm:
% 1. Determine representative width in longitunal direction of runup (uniform width)
% 2. For each runup ray determine for which Radius (R) water is encountered based on refelctivity or last point (A) -> store that R
% 3. Link time to each runup slice with the angular velocity info
% 4. Plot wet Radius R_wet (y-axis) v.s. time (x-axis) (2D): most probably an oscillation of R

% Detailed script for waterline algorithm:
% 1. Import 3D matrix
% 2. for loop through each ray
% 3. In for loop make if statement: wet or dry? 
% 4. If wet: store RTPA value in martix (wet_points) (rows = # rays, columns = RTPA) of that specific wet point in the ray 
% 5. Plot wet_points matrix (thick Red Line) along with your normal points 
% 6. Save image
% 7. repeat 1-7 but for different if statements threshold value for A as a sensitibity analysis

%% Summary:
%{ 
0. Supply [INPUT] parameters
1. Import dataset (RTPAXYZ_rays)
2. Use arc settings to obtain time information 
3. Noise reduction -> Manually or CloudCompare 
4. 3D to 2D matrix: RTPAXYZt or only RZt
5. Structred RZt 
6. Wave time-stack 
7. Wave Animation 
8. Export time-series data
-------------------TODO-----------------------
9. Wave Slider Plot [X]
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
runAnimation = false; % [INPUT] for running animation
saveOn2 = false; % [INPUT] boolean for saving gif's

% RUN ID, Date, Time, Resolution, angular range, scan duration, angular increment
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

%% 1. Import data
struct = load(['../Data/mat/RTPAXYZ_rays/',GPS_dataset(4:8),'/',map_dataset,'.mat']);
%struct = load(['../Data/mat/RTPAXYZ_rays/',GPS_dataset(4:8),'/',map_dataset,'_V2','.mat']);
RTPAXYZ_rays = struct.RTPAXYZ_rays;
clear struct;

addpath('../Functions/');
%% 2. Use arc settings to obtain time information
rays = size(RTPAXYZ_rays,3);
t_dur = str2double(MetaData{index,6}(1:2))*60+str2double(MetaData{index,6}(4:5)); % Total scan duration in [s]
t_delta = t_dur*MetaData{index,7}/MetaData{index,5}; % time between two rays
%t_delta = 0.1;

%% 3. Noise reduction:
%%{
% Also incl. gradient filtering later on maybe 
switch GPS_dataset
    case 'LW_DAY_3'
        filter = 1.3;
    case 'LW_DAY_4'
        filter = 2.0;
    case 'LW_DAY_5'
        filter = 2.0;
end

if count == 0
    % Filtering criterion: 
    A = RTPAXYZ_rays(:,7,:) > nanmedian(nanmedian(RTPAXYZ_rays(:,7,:))) + 5*nanstd(nanstd(RTPAXYZ_rays(:,7,:)));
    %A = RTPAXYZ_rays(:,7,:) > filter;
    B = [A, A, A, A, A, A, A]; % Apply filtering boolean for all columns
    RTPAXYZ_rays(B) = NaN; % Filter rows out 
    %RTPAXYZ_rays(:,7,:) = medfilt1(RTPAXYZ_rays(:,7,:),100,'omitnan');    
    count = -999;
end
%}

%% 4. 3D to 2D matrix: RTPAXYZt or only RZt 
RZt = zeros(size(RTPAXYZ_rays,1)*rays,3); % to hold data for Time-stack: R, Z and t
R = zeros(size(RTPAXYZ_rays,1)*rays,1);
Z = zeros(size(RTPAXYZ_rays,1)*rays,1);
t = zeros(size(RTPAXYZ_rays,1)*rays,1);

count_j = 1; % used to keep track of the row number 

for i=1:rays
    for j=1:size(RTPAXYZ_rays,1)
        if ~isnan(RTPAXYZ_rays(j,1,i)) && (RTPAXYZ_rays(j,1,i))~=0
            R(count_j,1) = RTPAXYZ_rays(j,1,i);
            Z(count_j,1) = RTPAXYZ_rays(j,7,i);
            t(count_j,1) = i;
            count_j = count_j + 1;
        end
    end
end
R(R == 0) = [];
Z(Z == 0) = [];
t(t == 0)= [];
RZt = [R,Z,t]; % Combine all columns in 2D matrix

%% 5. Structured RZt 
R_min = 20; 
R_max = 60;
R_delta = 0.25;

R_struc = R_min:R_delta:R_max; % Scattered R points will be interpolated on R_struc 
t_struc = unique(RZt(:,3)); % Corresponding t values of points
Z_struc = zeros(length(R_struc),length(t_struc));  % Corresponding Z values of points

count_t = 1; % used to keep track of ray number (and hence time, t)
count_rz = 1; % used to keep track of clustered points (RZ) of one ray
RZ_temp = []; % will hold RZ values of one ray temporary 


for i=1:size(RZt,1) % loop trough all (3D) points
    % cluster points of one ray:
    if RZt(i,3)==t_struc(count_t) 
        RZ_temp(count_rz,1) = RZt(i,1);
        RZ_temp(count_rz,2) = RZt(i,2);
        count_rz = count_rz + 1;
    % Increase ray number and interpolate and save previous Z values:
    else
        if nnz(RZ_temp~=0) >= 4 % needed since at least two points are needed for interpolation
           Z_struc(:,count_t) = (interp1(RZ_temp(:,1),RZ_temp(:,2),R_struc));
        else
           Z_struc(:,count_t) = ones(length(R_struc),1)*nan; % otherwise ray will contain just nan values
        end
        RZ_temp = []; % empty temprary RZ matrix for new ray
        count_rz = 1; % reset row number of RZ matrix for new ray
        count_t = count_t + 1; % increase ray number to next ray
        if RZt(i,3)==t_struc(count_t) % cluster points of this 1st next ray:
            RZ_temp(count_rz,1) = RZt(i,1);
            RZ_temp(count_rz,2) = RZt(i,2);
            count_rz = count_rz + 1;
        end
    end
end
t_0 = min(t_struc)*t_delta;
t_struc = (t_struc-min(t_struc))*t_delta;
Z_struc((Z_struc(:,:)==0)) = NaN;

%% 6.1 Wave Time-Stack:
t_min = 1; % [INPUT] start of time range in ray numbers
t_max = size(t_struc,1); % [INPUT] end of time range in ray numbers

Z_struc2 = Z_struc - min(min(Z_struc));

N =1000;
RGB = [0 150; 50 255; 200 250];
for n=1:N
    blue(n,:) = [RGB(1,1)+n*(RGB(1,2)-RGB(1,1))/N, RGB(2,1)+n*(RGB(2,2)-RGB(2,1))/N, RGB(3,1)+n*(RGB(3,2)-RGB(3,1))/N]/255;
end
colors = {gray,blue,parula,flag,winter,cool,jet}; % store colormaps

figure('units','normalized','outerposition',[0 0 1 1]);
[X,Y] = meshgrid(R_struc, t_struc);
colormap(colors{1});
contourf(Y(t_min:t_max,:),X(t_min:t_max,:),Z_struc(:,t_min:t_max)',10); hold on;
%contourf(Y(t_min:t_max,:),X(t_min:t_max,:),Z_struc2(:,t_min:t_max)',10); hold on;
h = colorbar;
set( h, 'YDir', 'reverse' );
shading flat;
%set(h, 'ylim', [1 2]);
title(['Wave time-stack (', map_dataset ,')']);
ylabel('Radius [m]');
xlabel('Time [s]');
%xlim([0, 100]); ylim([20,60]);
%xlim([0, 360]); ylim([40,50]);
h.Label.String = 'Elevation [m]';
set(gca,'fontsize',24);

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_RZt.png'],'-dpng');
    close all;
end

%% 6.2 Wave Snap-Shot:
% stores snap-shots:
snap = [786,42, 145, 205, 247, 372, 548, 623, 676, 925];  
snap2 = [50, 80, 160, 210];
snap3 = [50, 54, 58, 62, 64];

step = 6; % [INPUT] time-step in rays for plot 
a = round(11*3);
b = a + 4*step;
snap3 = a:step:b;

clear leg;
figure();

z_min = min(min(Z_struc(:,160:186)))-0.1;
z_max = max(max(Z_struc(:,160:186)))+0.3;

for i=1:length(snap3)
    plot(R_struc,Z_struc(:,snap3(i)),'-.','linewidth',3); hold on;
    if i == 1
        leg{i} = ['t_0 = ', num2str(round(t_struc(snap3(1)),1)) ,'s'];
    else
        leg{i} = ['t = t_0 + ', num2str(round((i-1)*t_delta*step,1)),'s'];
    end
end

xlim([R_min,R_max]); ylim([z_min,z_max]);
xlim([30,55]); ylim([0.2,1.6]);
set(gca,'fontsize',14);
title(['Wave development (', map_dataset,')']);
xlabel('Radius [m]');
grid;
ylabel('Water level elevation [m]');
legend(leg);

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_RZ.png'],'-dpng');
    close all;
end

%% 6.3 Wave Time-series:
r = 45;
r_index = find(R_struc==r); % [INPUT] index of radius of time-series
z_min = min(Z_struc(r_index,:))-0.1;
z_max = max(Z_struc(r_index,:))+0.1;

%Z_struc = detrend(Z_struc);
%A = detrend(Z_struc);

figure();
plot(t_struc,Z_struc(r_index,:),'linewidth',1.5);
%ylim([z_min,z_max]);
set(gca,'fontsize',14);
grid;
title(['Wave time-series at R = ', num2str(r),'m (', map_dataset,')']);
ylabel('Water level elevation [m]');
xlabel('Time [s]');

if saveOn
    print(['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'_Zt.png'],'-dpng');
    close all;
end

%% 7. Wave Animation:
if runAnimation
    filename = ['../Data/Images/',GPS_dataset(4:8),'/',map_dataset,'.gif'];
    
    f = figure();
    %ax = axes('Parent',f,'position',[0.15 0.25 0.75 0.65]);
    
    t_start = 1; % [INPUT]
    t_end = size(t_struc); % [INPUT]
    remove_ray = 0; % used for animation
    
    for i=t_start:t_end
        %g(i) = plot(RTPAXYZ_rays(:,1,i), RTPAXYZ_rays(:,7,i),'r-.','linewidth',1.5); hold on;
        %h(i) = plot(RZt(RZt(:,3)==i,1), RZt(RZt(:,3)==i,2), 'b-.','linewidth',1.5);
        h(i) = plot(R_struc, Z_struc(:,i), 'b-.','linewidth',1.5);
        xlim([R_min,R_max]);ylim([0.4,3]);
        xticks(linspace(R_min,R_max,5));
        title(['Animation: t=', num2str(round(i*t_delta,1)),'s (', map_dataset, ')']);
        xlabel('Radius [m]');
        ylabel('Elevation [m]');
        set(gca, 'fontsize', 14);
        
        pause(0.001);
        if saveOn2
            %Capture the plot as an image
            j=i-t_start+1;
            M(j) = getframe(f);
            im = frame2im(M(j));
            [imind,cm] = rgb2ind(im,256);
            %Write to GIF File
            
            if j==1
                imwrite(imind,cm,filename,'gif','Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append')
            end
        end
        if i > remove_ray
            delete(h(i-remove_ray));
            %delete(g(i-remove_ray));
        end
    end
    h(i) = plot(RTPAXYZ_rays(:,1,i), RTPAXYZ_rays(:,7,i),'b-.','linewidth',1.5); hold on;
    xlim([R_min,R_max]);ylim([0.4,3]);
    xticks(linspace(20,50,7));
    title(['Animation: t=', num2str(round(i*t_delta,1)),'s (', map_dataset, ')']);
    xlabel('Radius [m]');
    ylabel('Elevation [m]');
    set(gca, 'fontsize', 14);
end

%% 8. Export time-series data:
if saveOn
    save(['../Data/mat/RZt/',GPS_dataset(4:8),'/',map_dataset],'R_struc','Z_struc','t_struc','t_0');
end

%% 9. Wave Slider Plot:
%ray_num = 1;
%bgcolor = f.Color;
%b = uicontrol('Parent',f,'Style','slider','Position',[75,10,435,23],'value',ray_num, 'min',1, 'max',rays);
%bl1 = uicontrol('Parent',f,'Style','text','Position',[50,12,23,23],'String','0','BackgroundColor',bgcolor);
%bl2 = uicontrol('Parent',f,'Style','text','Position',[510,12,28,22],'String',num2str(rays),'BackgroundColor',bgcolor);
%bl3 = uicontrol('Parent',f,'Style','text','Position',[245,35,100,20],'String','Ray number','BackgroundColor',bgcolor);
%b.Callback = @(es,ed) updateSystem(h,ray_num);
