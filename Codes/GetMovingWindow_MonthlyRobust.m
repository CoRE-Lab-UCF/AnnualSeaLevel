function [RTotalPreMSL] = GetMovingWindow_MonthlyRobust(MSL,TIME,Window)
% MSL : mean sea level time series, vector (n*1)
% TIME, vector(n*1)
Window = 5; % window: scalar (years) , like:  5 
% output: a structure variable
% we do analysis in each moving window, and the window is shifted by one month each time step.
% for example, we set Window=5 year. then in each fitting analysis (like 2000-2004),
% the Amplitude and phase is regarded as the results at  
% the central month of the time series, namely, the 30th month (2002-06).
 
m = length(MSL);
CONT = 0;
PreMSL = [];
PreSASSA = [];
for j = 1:m-Window*12+1
    
    L = MSL(j:j+Window*12-1);   % build matrix
  
%%  Time datum based on real time       
    t = [TIME(j):1/12:TIME(j+Window*12-1)]';
    A(1:Window*12,1) = ones(Window*12,1);
    A(1:Window*12,2) = t/Window*12;
    A(1:Window*12,3) = sin(2*pi/0.5*t);
    A(1:Window*12,4) = cos(2*pi/0.5*t);
    A(1:Window*12,5) = sin(2*pi/1*t);
    A(1:Window*12,6) = cos(2*pi/1*t);
    
    AA = A;
    BOOL = isnan(L);
    if sum(BOOL)>=12*2  % set threshold stop fitting (if missing data is too much) 
        
        if j == 1  % the beginning of the fitting
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
            PreMSL(j,3:10) = NaN;
        end
        AmpPha(j,1) = TIME(Window*6+j-1,1);
        AmpPha(j,2:5) = NaN;                
    else  
        A(BOOL,:) = [];
        L(BOOL) = [];        
%         X = inv(A'*A)*A'*L;  % least squares analysis
%         V = A*X-L;
%         sigma = sqrt(V'*V/(Window*12-sum(BOOL)-6));  % std
%         DX = sqrt(diag(sigma^2*inv(A'*A))); % std of each parameter
%         t = abs(X)./DX; % significance test (one sigma level)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%   Robust estimation        
        [X STAT]= robustfit(A(:,2:6),L,'huber');  
        % 1) 'huber' is one kind of weight function in robust estimation to
        %    weaken the effect of outliers on results; You can choose other functions
        %    or use default function. Regarding default function: 
        %             [X STAT]= robustfit(A(:,2:6),L); 
        %    For other weight functions:
        %             [X STAT]= robustfit(A(:,2:6),L,'bisquare');               
        % 2) STAT is a struct which saves statistical results of robust
        %    estimation; for more detail, please refer to 'robustfit'
        %    introduction.
        % 3) You may encounter "Warning: Iteration limit reached.". This means
        %    the robust estimation is not converged. Try to open function "statrobustfit.m" 
        %    and go to line 75 where is coded with "iterlim = 50;". Change 50
        %    to 100, 1000 or 10000, and run robust estimation again to see the results. 
     
        sigma = STAT.robust_s;
        DX = STAT.se;
        t = abs(X./STAT.se);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
       
        TEMP = AA(:,1:2)*X(1:2);  % reconstruct trend
        PreMSL(j,1) = TIME(Window*6+j-1,1);% + TIME(Window*6+j-1,2)/12-1/12;
        PreMSL(j,2) = TEMP(Window*6);        
        TEMP = AA(:,3:4)*X(3:4); % reconstruct semi-annual cycle
        PreMSL(j,3) = TEMP(Window*6);               
        TEMP = AA(:,5:6)*X(5:6); % reconstruct annual cycle
        PreMSL(j,4) = TEMP(Window*6);
        PreMSL(j,5) = norm(X(3:4));   %  semi-annual amplitude
        PreMSL(j,6) = norm(X(5:6));   %  annual amplitude

        DXX = sigma^2*inv(A'*A);
        TEMP1 = [X(3)/sqrt(X(3)^2 + X(4)^2), X(4)/sqrt(X(3)^2 + X(4)^2)];
        TEMP2 = [X(5)/sqrt(X(5)^2 + X(6)^2), X(6)/sqrt(X(5)^2 + X(6)^2)];                                   
              
        PreMSL(j,7) = sqrt(TEMP1*DXX(3:4,3:4)*TEMP1');  % std of semi-annual amplitude
        PreMSL(j,8) = sqrt(TEMP2*DXX(5:6,5:6)*TEMP2');  % std of annual amplitude
        PreMSL(j,9) = sqrt(AA(Window*6,3:4)*DXX(3:4,3:4)*AA(Window*6,3:4)'); % std of semi-annual cycle
        PreMSL(j,10) = sqrt(AA(Window*6,5:6)*DXX(5:6,5:6)*AA(Window*6,5:6)');  % std of annual cycle
        PreErr = [PreMSL(j,9) PreMSL(j,10)];
        if j==1
            Beginning(1:Window*6-1,1) = TIME(1:Window*6-1,1);
            TEMP = AA(:,1:2)*X(1:2);
            Beginning(1:Window*6-1,2) = TEMP(1:Window*6-1,1);
            TEMP = AA(:,3:4)*X(3:4);
            Beginning(1:Window*6-1,3) = TEMP(1:Window*6-1,1);
            TEMP = AA(:,5:6)*X(5:6);
            Beginning(1:Window*6-1,4) = TEMP(1:Window*6-1,1);
        elseif j==m-Window*12+1
            Ending(1:Window*6,1) = TIME(end-Window*6+1:end,1);
            TEMP = AA(:,1:2)*X(1:2);
            Ending(1:Window*6,2) = TEMP(Window*6+1:Window*6*2,1);
            TEMP = AA(:,3:4)*X(3:4);
            Ending(1:Window*6,3) = TEMP(Window*6+1:Window*6*2,1);
            TEMP = AA(:,5:6)*X(5:6);
            Ending(1:Window*6,4) = TEMP(Window*6+1:Window*6*2,1);     
        end        
        %%
        AmpPha(j,1) = TIME(Window*6+j-1,1);
        AmpPha(j,2) = norm(X(3:4));
        AmpPha(j,3) = norm(X(5:6));
        [ang] = GetPhase2(X(3),X(4));
        AmpPha(j,4) = ang/pi*180;       % semi-annual phase
        [ang] = GetPhase2(X(5),X(6));     
        AmpPha(j,5) = ang/pi*180;       % annual phase         
    end 
end
RTotalPreMSL.Monthly = [Beginning,zeros(Window*6-1,6);PreMSL;Ending,zeros(length(Ending),6)];

% RTotalPreMSL.Monthly. Column 1: Time (years), 
%                       Column 2: Secular Trend ,  
%                       Column 3: semi-annual cycle,  
%                       Column 4: annual cycle,   
%                       Column 5: semi-annual amplitude,   
%                       Column 6: annual amplitude,   
%                       Column 7: standard error of  semi-annual  amplitude ,  
%                       Column 8: standard error of  annual amplitude,   
%                       Column 9: standard error of  semi-annual cycle,  
%                       Column 10: standard error of  annual cycle.

RTotalPreMSL.MSL = MSL;
RTotalPreMSL.TIME =TIME;
RTotalPreMSL.AmpPha = AmpPha;
% AmpPha:   Column 1: Time (years)
%           Column 2: semi-annual amp
%                  3: annual amplitude
%                  4: semi-annual phase
%                  5: annual phase





