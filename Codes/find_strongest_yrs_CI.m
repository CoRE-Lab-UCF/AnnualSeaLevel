function [yrs5_nino34,yrs5_LAnina34,no_overlap] = find_strongest_yrs_CI(names,YRSelnino34,first12, last12, n34_Ann_Avg_nan2, data10)

%% Find Strongest Years of Climate Index

% Pre-allocate index for beginning and end of time
beg = zeros([length(names) 1]); % what time we begin to check
fin = zeros([length(names) 1]); % what time we stop checking

no_overlap = [];

% For each tide gauge, take the Annual Maximum during the alotted time period
for tg = 1:length(names) % for each tg

    yrsTG.(names{tg}) = unique(fix(data10.(names{tg})(:,1))); % list of yrs for each TG
    
    % Calculate the length of overlapping years
    % Find index in Climate TS 
    [iCI{tg}, loc_CI{tg}] = ismember(yrsTG.(names{tg}), YRSelnino34.(names{tg})); % find index in NAO TS

 % If TG and Climate Index overlap for less than 30 years
    if sum(iCI{tg}) < 30  

    % If tg and CI do not overlap 
    %if last12(tg)<= YRSelnino34.(names{tg})(1)

        no_overlap = cat(1,no_overlap,tg);
        test5years(tg, 1:5) = NaN(1,5);
        lanina5yrs(tg, 1:5) = NaN(1,5);

    else


        if first12(tg) <= YRSelnino34.(names{tg})(1); % If tg starts earlier than climate index record,
            beg(tg) = 1;              % Begin at 1st yr of climate index record
        else        % If tg starts after the 1st yr of climate index record
            beg(tg) = find(YRSelnino34.(names{tg}) == first12(tg)); % begin at the 1st yr of tg record
        end
        if last12(tg) == 2019; % the originally downloaded TG data ends in 2019
            fin(tg) = find(YRSelnino34.(names{tg}) == 2019);   % Finish checking until this time period
        else                   % If the TG MSL record ends before 2019...
            fin(tg) = find(YRSelnino34.(names{tg}) == last12(tg)); % Finish scanning climate index until the last yr of TG MSL record
        end

        % For each tide gauge, take the top 5 Annual Maximum of the Climate Index
        % during the alotted time period (time of TG MSL record)
        [nPos(tg,1:5), iNp(tg,1:5)] = maxk(n34_Ann_Avg_nan2.(names{tg})(beg(tg):fin(tg)),5); % This index is based off the shortened TS...
        % Align the index properly (iNp refers to the shortened section of the TS)
        yrs_tg = YRSelnino34.(names{tg})(beg(tg):fin(tg));
        test5years(tg, 1:5) = yrs_tg(iNp(tg,:));


        % Also find La Niña years (by taking minimum)
        [nNeg(tg,1:5), iNneg(tg,1:5)] = mink(n34_Ann_Avg_nan2.(names{tg})(beg(tg):fin(tg)),5);
        % Align index properly
        lanina5yrs(tg,1:5) = yrs_tg(iNneg(tg,:));

    end

        % Sort the years in ascending order
        yrs5_nino34 = sort(test5years, 2); % sort along 2nd dimension (rows)
        yrs5_LAnina34 = sort(lanina5yrs,2);


end