%% This is the main program to run the seasonal sea level cycle analysis.

clc  % clear command window
clear % clear variables

%% Run results for all tg's in the threshold (3 amp diff & .8 RMSE bt models)

load TGdata38.mat
name = {TGdata38.name}';
% Remove special characters(_)from field names
for j = 1:length(name)  
    TEMP = name{j};
    BOOL = find(TEMP==' ');
    TEMP(BOOL) = '_';
    name{j} = TEMP;
end
valid_names = matlab.lang.makeValidName(name); % default underscore, but doesn't underscore spaces.

tic  
for j = 1:length(TGdata38)   
    MSL = TGdata38(j).height;        
    % Convert Time from current Matlab format to PSMSL decimal format 
    jd = TGdata38(j).time;         %(time in Julian Day)
    time = datevec(jd);
    TIME = time(:,1) + (time(:,2)/12) - 1/24; 
%        year     +       month    - day(24 hrs)
    Window = 5;
    [RTotalPreMSL] = GetMovingWindow_MonthlyRobust(MSL,TIME,Window); % Using Robust Fit 12 mins)
 
    RobustFitResults.(valid_names{j}) = RTotalPreMSL;               % save shorter list of the above (code using atan(a/b), c=a, ang=ang)
end 
toc
