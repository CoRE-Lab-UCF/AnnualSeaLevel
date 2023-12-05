%% CORRECTED Climate Index Analysis

clear; close all; clc

% For each climate index, check the 5 strongest years in both the positive and negative phase.
% Apply the Harmonic analysis to the MSL during these years
% Compare the output: amplitude value in (+) and in (-)
% BoxPlot the difference between these with standard errors

% Not all TGs have data during the strongest years of a climate index...
% Align time of Climate Index w TG time frame

cd('C:\Users\am612589\OneDrive - University of Central Florida\DataCode_Inc\Supplementary')
addpath 'C:\Users\am612589\OneDrive - University of Central Florida\DataCode_Inc\Supplementary\IndicesPSL_NOAA';
addpath('C:\Users\am612589\OneDrive - University of Central Florida\DataCode_Inc\Supplementary\functions')
addpath('C:\Users\am612589\OneDrive - University of Central Florida\DataCode_Inc\cbarrow_v1.1\cbarrow')

%% Align Time to TG MSL record
load merged_data;
load ResMerged;
%load ENSO34ts; % Time vector for Harmonic Analysis

names = fieldnames(ResMerged);

%% Load Climate Index: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NAO
% load NAOdata
% load NAO_Ann_Avg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AMO
% load AMOdata
% load AMO_Ann_Avg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDO
% load PDOdata
% load PDO_Ann_Avg

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% El Niño 1.2
%n12data = table2array(readtable('Nino 1 2 anomaly 1987-2010 mean removed.txt'));
% 
% Turn monthly data into annual average
%n12_Ann_Avg = nanmean(n12monthly,1)
% load n12_Ann_Avg
% load n12yrs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% El Niño 3 Anomaly
% Turn monthly data into annual average
%n3_Ann_Avg = nanmean(ENSO3monthly,1)
load n3_Ann_Avg
load ENSO3yrs

% Change variable names
Ann_Avg_IC = n3_Ann_Avg
yrs = ENSO3yrs; % YRS of climate index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preprocessing
[first12,last12,YRSCI,Ann_Avg_nan2,data10] = pre_process_index(ResMerged,yrs,Ann_Avg_IC);

%% Find Strongest Years of Climate Index
[yrs5_nino34,yrs5_LAnina34,no_overlap] = find_strongest_yrs_CI(names,YRSCI,first12, last12, Ann_Avg_nan2, data10);

yrs5_nino34(isnan(yrs5_nino34(:,1)),:)     = [];
yrs5_LAnina34(isnan(yrs5_LAnina34(:,1)),:) = [];

data10 = rmfield(data10,names(no_overlap));

%% Subset MSL during Climate Index 5 Strongest Yrs
% Save for each TG
MSL5_CI = [];
names = fieldnames(data10);

for tg = 1: size(yrs5_nino34,1) % through tide gauges

    MSL_y = [];

    % top 5 Climate Index years for each tg 
    year_tg = yrs5_nino34(tg,:);

    a = fix(data10.(names{tg})(:,1)); % all years with >10 data months/yr

    for ii = 1: length(year_tg)      % For each of the 5 strongest years 

        fa = find(a== year_tg(ii));
        MSL_y = cat(1,MSL_y,[data10.(names{tg})(fa,1), data10.(names{tg})(fa,2)]); % concatenate the results as a timeseries

    end

    MSL5_CI.(names{tg}) = MSL_y;

end


%% Run Harmonic Analysis El Niño 3.4
% Create Time Vector
% The harmonic needs a time vector. The 5 strongest years are not truly chronological so just
% use some 5-year period for the timestep to reference
% Let's use the first 5 years

TIME = ENSO34ts(1:5*12,1);
%%% El Niño 3.4
for i = 1:length(names); 
    x = MSL5_CI.(names{i})(:,2); % MSL Values during the 5 strongest yrs of El Niño 3.4  
    
    Window = 5;
    [res_nino34] = GetMovingWindow_MonthlyRobust(x,TIME,Window);  % Using Robust Fit 
    res_MSL_CI.(names{i}) = res_nino34;                     % save Robust Fit results to a struct  
end

%% Negative Phase: La Niña
clear year_tg
clear fa
MSL5_LaNina34 = [];

for tg = 1: size(yrs5_LAnina34,1) % through tide gauges

    MSL_ynina = [];

    % top 5 Climate Index years for each tg 
    year_tg = yrs5_LAnina34(tg,:);

    a = fix(data10.(names{tg})(:,1)); % all years with >10 data months/yr

    for ii = 1: length(year_tg);      % For each of the 5 strongest years 

        fa = find(a== year_tg(ii));
        MSL_ynina = cat(1,MSL_ynina,[data10.(names{tg})(fa,1), data10.(names{tg})(fa,2)]); % concatenate the results as a timeseries

    end

    MSL5_LaNina34.(names{tg}) = MSL_ynina;

end

%% Run Harmonic Analysis: La Niña
% Time vector will be the same (bc it doesn't matter)
for i = 1:length(names); 
    x = MSL5_LaNina34.(names{i})(:,2); % MSL Values during the 5 strongest yrs of El Niño 3.4  
    
    Window = 5;
    [res_nina34] = GetMovingWindow_MonthlyRobust(x,TIME,Window);  % Using Robust Fit 
    res_MSL_LaNina34.(names{i}) = res_nina34;                     % save Robust Fit results to a struct  
end

%% Compute diff b/t + & - phase of El Nino
for k = 1:length(names)
    n34_diff(k) =  res_MSL_CI.(names{k}).AmpPha(3) - res_MSL_LaNina34.(names{k}).AmpPha(3);
end

%% Map the difference
load lon_merged
load lat_merged

lon_merged(no_overlap) = [];
lat_merged(no_overlap) = [];

% hh= figure;
% 
% set(hh,'units','inches','Position', [0.05    0.02    7.6   3.5 ],'InvertHardCopy','off',...
%     'resize','off','PaperPositionMode','auto','PaperType','A0','visible','on',...
%     'color','w');
% 
% m_proj('miller','long', [-179.8667 176.6811],'lat',[-66.7462 81.1167] );  % projection
% m_gshhs('ir','patch',[0.90 .90 .90],'edgecolor','w');hold on   % show Lakes
% % m_scatter(LON,LAT,Y,S,C)
% m_scatter(lon_merged, lat_merged, 10, abs(n34_diff),  'filled') ; hold on 
% m_grid('tickdir','in','linest','none','FontName','Times','FontSize',12,'xticklabels',[],'yticklabels',[]);       
% 
% title('Abs Diff b/t Amp During El Niño & La Niña 5 Strongest Yrs')
% colormap('viridis')
% cb = colorbar();
% 
% caxis( [min(abs(n34_diff)) 6])
% cbarrow('up')
% cb.Label.String = '(mm)'

%% Location example
% Boxplot Alex's version

% Somewhere in the Pacific where ENSO impacts.. 
% Honolulu
% HONOLULU
xx = [res_MSL_CI.HONOLULU.AmpPha(3) ;
    res_MSL_CI.HONOLULU.AmpPha(3) + res_MSL_CI.HONOLULU.Monthly(end,8);
    res_MSL_CI.HONOLULU.AmpPha(3) - res_MSL_CI.HONOLULU.Monthly(end,8)];

yy = [res_MSL_LaNina34.HONOLULU.AmpPha(3) ;
    res_MSL_LaNina34.HONOLULU.AmpPha(3) + res_MSL_LaNina34.HONOLULU.Monthly(end,8);
    res_MSL_LaNina34.HONOLULU.AmpPha(3) - res_MSL_LaNina34.HONOLULU.Monthly(end,8)];

figure; hold all;
boxplot([xx,yy],'Labels',{'Postive phase','Negative phase'})
title('Honolulu')
grid minor
ylabel('Mean annual amplitude (mm)')

%% Significance of amp diff
% Compute the significance of the amplitude difference during the Positive Phase vs. the Negative Phase
% Check if the error bars overlap

names = fieldnames(data10);
amp5_nino34 = zeros([3, length(names)]); % Pre-allocate for el nino amp, & SE
amp5_nina34 = nan([3,length(names)]); 

for tg = 1:length(names)
    amp5_nino34(1:3 , tg) = [res_MSL_CI.(names{tg}).AmpPha(3); 
                             res_MSL_CI.(names{tg}).AmpPha(3) + res_MSL_CI.(names{tg}).Monthly(end,8);
                             res_MSL_CI.(names{tg}).AmpPha(3) - res_MSL_CI.(names{tg}).Monthly(end,8);] ;

    amp5_nina34(1:3 , tg) = [res_MSL_LaNina34.(names{tg}).AmpPha(3); 
                             res_MSL_LaNina34.(names{tg}).AmpPha(3) + res_MSL_LaNina34.(names{tg}).Monthly(end,8);
                             res_MSL_LaNina34.(names{tg}).AmpPha(3) - res_MSL_LaNina34.(names{tg}).Monthly(end,8);] ;

    if amp5_nino34(1,tg) > amp5_nina34(1, tg) % If El Niño amp > La Niña amp
        sig(tg) = amp5_nino34(3,tg) > amp5_nina34(2,tg); % Compare the bottom errorbar of El Niño to the top errorbar of La Niña

    else 
        sig(tg) = amp5_nina34(3,tg) > amp5_nino34(2,tg); 
        
    end

end

%% Map significant values

hh= figure;

set(hh,'units','inches','Position', [0.05    0.02    7.6   3.5 ],'InvertHardCopy','off',...
    'resize','off','PaperPositionMode','auto','PaperType','A0','visible','on',...
    'color','w');

m_proj('miller','long', [-179.8667 176.6811],'lat',[-66.7462 81.1167] );  % projection
m_gshhs('ir','patch',[0.90 .90 .90],'edgecolor','w');hold on   % show Lakes
% m_scatter(LON,LAT,Y,S,C)
m_scatter(lon_merged(find(sig ==1)), lat_merged(find(sig ==1)), 12, n34_diff(find(sig==1)),  'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]) ; hold on 
m_grid('tickdir','in','linest','none','FontName','Times','FontSize',12,'xticklabels',[],'yticklabels',[]);       

title('Sig Amp Diffs During 5 Strong Yrs of Pos Phase & Neg Phase', 'Fontsize', 13.4)
cb = colorbar()
caxis([-90 90])
cb.YTick = -90:30:90

cb.Label.String = 'mm';

% we want a blue-red diverging colormap, but without white in the center
whatt = colormap(m_colmap('diverging'))
%Find the color white range (0.9 0.9 0.9) and remove it
whatt(114:142, :) = [];
colormap(whatt)

cbarrow



