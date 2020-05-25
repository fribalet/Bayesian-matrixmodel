function [year_day]=find_yearday(dates)
    
    %solution posted on
    %http://www.mathworks.com/matlabcentral/answers/1756-converting-a-date-string-to-day-of-year
    
    dv  = datevec(dates);  % gives an N x 6 array
    dv  = dv(:, 1:3);   % only need first three pieces of time
    dv2 = dv;
    dv2(:, 2:3) = 0;    % set day and month to 0, corresponds to Dec 31st...
    year_day = datenum(dv) - datenum(dv2);

end
