%% Drop yrs of MSL TG data with <= 10 months

function [first12,last12,YRSelnino34,n34_Ann_Avg_nan2,data10] = pre_process_index(ResMerged,yrs,n34_Ann_Avg)

names = fieldnames(ResMerged);

% copy MSL data to drop years (& be able to access full data later)
for i = 1:length(names)

    data10.(names{i}) = [ResMerged.(names{i}).TIME, ResMerged.(names{i}).MSL]; % copy MSL data to drop years (& be able to access full data later)

    MSL10.(names{i}) = ResMerged.(names{i}).MSL; % simplify name of MSL variable

    years{i} = fix(ResMerged.(names{i}).TIME); % Years of TG data 

    % Find where there are more than 2 NaN months per year
    u_yrs{i} = unique(years{i}); % For each tg, List the unique years
    for u = 1:numel(u_yrs{i})
        iyr = u_yrs{i}(u); % For each yr
        nan_count_per_year{i}(u) = sum(isnan(MSL10.(names{i})  (years{i} == u_yrs{i}(u)) )) ;      
    end

    fnan{i} = find(nan_count_per_year{i} > 2); % Find where num of NaNs exceed our threshold
    
    % change positions to years
    yrsXX{i} = u_yrs{i}(fnan{i});
    % Pre-allocate space for indexing
    fn_yrs = []; % index for the years we need to drop in TG data
    iDCyrs = []; % index in Climate index for the yrs we need to drop (El Niño)
    
    for n = 1:length(fnan{i}) 
        % In the TG data, find the index of years we need to drop
        temp{i} = find(years{i} == yrsXX{i}(n)); % Find index for each year that has more than 2 NaNs
        fn_yrs = cat(1,fn_yrs,temp{i});          % Concatenate the indices for each tg and each yr 

        % In the climate index, find the index for the same (corresponding) yrs: El Niño
        Ctemp{i} = find(yrs == yrsXX{i}(n));
        iDCyrs = cat(1, iDCyrs, Ctemp{i});

    end

    % Drop the years
             data10.(names{i})(fn_yrs, :) = []; % Now this is the data (TIME,MSL) after removing the years w <10 months
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% El Niño         
             YRSelnino34.(names{i}) = yrs; % Create a structure for the yrs of the Climate Index associated with each TG
             YRSelnino34.(names{i})(iDCyrs) = []; % Drop the necessary years
             n34_Ann_Avg_nan2.(names{i}) = n34_Ann_Avg; % Create a structure for the Climate Index values during the years w >= 10 months of data
             n34_Ann_Avg_nan2.(names{i})(iDCyrs) = [];  % Drop the same years we dropped 

    % Make sure the first/last years we use have the full 12 months
    u12{i} = unique(fix(data10.(names{i})(1:12,1))); % Check if first 12 months are in the same year
    L(i) = length(u12{i});           % # of yrs displayed in first 12 months
  
    if L(i) == 2                     % If the first 12 months occurs b/t 2 years... 
        first12(i) = u12{i}(2);      % Consider the 2nd year (with full months) to be the first where we scan for the climate Index 
    else
        first12(i) = u12{i};
    end

    ulast{i} = unique(fix(data10.(names{i})(end-12+1:end,1)));
    Llast(i) = length(ulast{i});
    
    last12(i) = ulast{i}(1);   % Consider the yr with full months to be the last year we check          
end


end
