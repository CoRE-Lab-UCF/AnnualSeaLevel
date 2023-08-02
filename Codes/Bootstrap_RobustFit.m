function [bsResBox] = Bootstrap_RobustFit(MSL,TIME,Window)
% MSL : mean sea level time series, vector (n*1)
% TIME, vector(n*1)

Window = 5; % window: scalar (years), like: 5 
% output: a structure variable

% we do analysis in each moving window, and the window is shifted by one month each time step.
% for example, we set Window=5 year. then in each fitting analysis (like 2000-2004),
% the Amplitude and phase is regarded as the results at  
% the central month of the time series, namely, the 30th month (2002-06).

%% Sensitivity Test
% Leave-1-out Bootstrapping: for each 5-yr window we drop one year and fit the model
% (then drop the next year and fit the model, etc.) 
% this gives for each window 5 different estimates of the amplitudes and phases 
% we can look at the spread (range b/t upper & lower boundary) to see how individual years 
% influence the results. 

m = length(MSL);
CONT = 0;
PreMSL = [];
PreSASSA = [];
%RngMSL = zeros(j,12);  % Pre-allocate 
x_bootstrap = [];
y_bootstrap = [];
for j = 1:m-Window*12+1
    
    L = MSL(j:j+Window*12-1);   % build matrix
   
%%  Time datum based on real time       
    t = [TIME(j):1/12:TIME(j+Window*12-1)]';
    A(1:Window*12,1) = ones(Window*12,1);
    A(1:Window*12,2) = t/Window*12;
    A(1:Window*12,3) = sin(2*pi/0.5*t);   % semi-annual cycle
    A(1:Window*12,4) = cos(2*pi/0.5*t);   %             amplitude
    A(1:Window*12,5) = sin(2*pi/1*t);     % annual cycle
    A(1:Window*12,6) = cos(2*pi/1*t);     %             amplitude
    
    AA = A;
    BOOL = isnan(L);
    
    % 6 parameters for our regression model = 1(slope)+1(intercept) + 2(annual) + 2(semi-annual)

    if sum(BOOL)>=12*2  % set threshold stop fitting (if missing data is too much) 
        
        if j == 1  % the beginning of the fitting
            %Beginning(1:Window*6,1) = TIME(1:Window*6,1);
            %Beginning(1:Window*6,2:4) = NaN*ones(Window*6,3);
            Beginning(1:Window*6-1,1) = TIME(1:Window*6-1,1);
            Beginning(1:Window*6-1,2:4) = NaN*ones(Window*6-1,3);
            X = [1:6]';PreErr = [1 1];
            PreMSL(j,1) = TIME(Window*6+j-1);
            PreMSL(j,2) = NaN; 
            
        elseif j==m-Window*12+1   % the end of the fitting 
            Ending(1:Window*6+1,1) = TIME(end-Window*6:end,1);
            Ending(1:Window*6+1,2:4) = NaN*ones(Window*6+1,3);
            X = [1:6]';PreErr = [1 1];           
        else 
            X = [1:6]';PreErr = [1 1];       
            PreMSL(j,1) = TIME(Window*6+j-1);
            PreMSL(j,2) = NaN; 
            PreMSL(j,3:12) = NaN;

        end
        AmpPha(j,1) = TIME(Window*6+j-1,1);
        AmpPha(j,2:5) = NaN;
  
              
    else  
%         A(BOOL,:) = [];  % empty out the NaNs from A (x: time)  
%         L(BOOL) = [];    % remove NaNs from L (y: MSL)
%         
%         %% Build internal loop to have 5 estimates for each 5-yr. window %(j:j+Window*12-1)
%         for ij = j:j+Window*12-1  % each moving 5-yr window
%             if(ij~=j+12-1)        % reject one year of data
%                 x_bootstrap = [x_bootstrap A(ij)];    % store results of the model excluding 1 yr each time  
%                 y_bootstrap = [y_bootstrap L(ij)];
%             end
%         end

        % 5-year moving window means we have 60 (5*12=60) values in our window.
        % In my opinion, one-year bootstrap means we remove 12 values from the window and
        % run analysis.
        
        PreTemp = [];AmpPhaTemp = []; % set up Matrix to save results
        for k = 1:Window % we have five years data
            A0 = A;
            L0 = L;
            A0(12*(k-1)+1:12*k,:) = []; % remove 1 year data from the window
            L0(12*(k-1)+1:12*k,:) = [];
            BOOL = isnan(L0); % remove NaN
            L0(BOOL) = [];
            A0(BOOL,:) = [];
            [PreMSL,AmpPha] = RobFit(A0,L0,BOOL,A,TIME);  % RobFit is function for RobustFit estimation.    
            PreTemp = [PreTemp;PreMSL]; % save results
            AmpPhaTemp = [AmpPhaTemp;AmpPha];
        end
        
        A0 = A;  % caculate full 5-year window (like original code)
        L0 = L;
        BOOL = isnan(L0);
        L0(BOOL) = [];
        A0(BOOL,:) = [];

        [PreMSL,AmpPha] = RobFit(A0,L0,BOOL,A,TIME);    
        PreTemp = [PreTemp;PreMSL];  % save results for comparison
        AmpPhaTemp = [AmpPhaTemp;AmpPha];
        % Include the range of the bootstrap results to assess the spread 
        RngMSL = range(PreTemp(1:5,:));
        RngPha = range(AmpPhaTemp(1:5,:));
        
        PreTemp = [PreTemp;RngMSL];
        AmpPhaTemp = [AmpPhaTemp;RngPha];
        
        %MON = round( (TIME(j) - fix(TIME(j))) *12+1/24);
        MON = round( (TIME(j+ (Window*12/2 -1)) - fix(TIME(j+ (Window*12/2 -1)))) *12+1/24);
        STR = ['bs_',num2str(fix(TIME(j+ (Window*12/2 -1)))),'_',num2str(MON)]; % save as Struct
        bsResBox.(STR).PreMSL = PreTemp;
        bsResBox.(STR).AmpPha = AmpPhaTemp; 
        bsResBox.Range = RngMSL;
        %SizeRng = zeros(length(fields(ResBox))-1,12);
        %SizeRng(1:length(fields(ResBox))-1,1:12) = RngMSL;
        %ResBox.Range = [SizeRng]; 
        
        %% PreTemp. Size: 6*12. First 1-5 lines saves the results of
        % bootstrap. The 6th line saves the results of 5-year window
        % The same as AmpPhaTemp
        % I (Amanda) added the 7th row: range of the results throughout the 5-yr
         
        %= nanmax(RngMSL)
        %bsMax_Rng_AC(j,1) = nanmax(BootstrapRes.(STR).PreMSL(7,4)/10);  % range (row 7) of Annual Cycle (cm) col. 4
    end
 
end

%TotalPreMSL.Monthly = [Beginning,zeros(Window*6-1,6);PreMSL;Ending,zeros(Window*6,6)];
% TotalPreMSL.Monthly = [Beginning,zeros(Window*6-1,8);PreMSL;Ending,zeros(length(Ending),8)];

% TotalPreMSL.Monthly.  Column 1: Time (years), 
%                       Column 2: Secular Trend ,  
%                       Column 3: semi-annual cycle,  
%                       Column 4: annual cycle,   
%                       Column 5: semi-annual amplitude,   
%                       Column 6: annual amplitude,   
%                       Column 7: standard error of  semi-annual  amplitude ,  
%                       Column 8: standard error of  annual amplitude,   
%                       Column 9: standard error of  semi-annual cycle,  
%                       Column 10: standard error of  annual cycle.
%                       Column 11: std err of semi-annual phase
%                       Column 12: std err of annual phase

%% TotalPreMSL.Hourly = PreSASSA;
% TotalPreMSL.MSL = MSL;
% TotalPreMSL.TIME =TIME;
% TotalPreMSL.AmpPha = AmpPha;
% AmpPha:   Column 1: Time (years)
%           Column 2: semi-annual amp
%                  3: annual amplitude
%                  4: semi-annual phase
%                  5: annual phase 