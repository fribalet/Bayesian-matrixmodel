% make average solar dawn/dusk for a year:

%Can upload from preformed Solar.mat files, but if don't have that, can
%process raw data for both dawn and dusk hours:

temp_dawn=nan(366,11);
temp_dusk=nan(366,11);

envpath='/Volumes/Lab_data/MVCO/EnvironmentalData/';

for year2do=[2003:2009 2011:2015]
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    %if have premade Solar.mat files:
    %     eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year) '/model_code/solar.mat'])
    %     yrdays=find_yearday(dawn(:,1));
    %     temp_dawn(yrdays,year-2002)=dawn(:,2);
    
    %Otherwise, use the raw files:
    
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99'',''ascii'');'])
    
    yd_met=X(:,2);
    Year = X(:,1);
    Day = X(:,4);
    Hour = X(:,5);
    Solar=X(:,14);  %Solar_campmt_median
    
    clear X
    
    t = find(diff(yd_met) < 0);
    yd_met(t+1) = NaN;
    Solar(t+1) = NaN;
    
    date_met = yd_met + datenum(['1-0-' num2str(year2do)]);  %UTC, matlab date (change 11/4/03, also in setup_days.m)
    dawnlevel = 5;
    
    ind = find(~isnan(Solar) & ~isnan(yd_met));
    date_met = date_met(ind);
    yd_met = yd_met(ind);
    Solar = Solar(ind);
    Hour = Hour(ind);
    
    unqday = unique(floor(date_met - 5/24)); %approx number of local days
    
    unqday = unqday(~isnan(unqday));
    dawnhr=nan(size(unqday));
    duskhr=nan(size(unqday));
    for count = 1:length(unqday),
        ind = find(floor(date_met - 5/24) == unqday(count));
        ind2 = find(Solar(ind) > dawnlevel);
        if ~isempty(ind2),
            %        dawnhr(count) = Hour(ind(ind2(1))) - 1 - 4;  %local time from May onwards through summer, 1 h before dawnlevel?
            %        duskhr(count) = Hour(ind(ind2(end))) - 4 + 1;  %next hour after end of light
            dawnhr(count) = Hour(ind(ind2(1))) - 1;  %1 h before dawnlevel?
            duskhr(count) = Hour(ind(ind2(end))) + 1;  %next hour after end of light
            if dawnhr(count) == -1, keyboard, end;
        else
            dawnhr(count) = NaN;
            duskhr(count) = NaN;
        end;
    end;
    
    if year2do==2003;
        dawnhr=dawnhr(2:end)';
        duskhr=duskhr(2:end)';
        unqday=unqday(2:end)';
    end
    yrdays=find_yearday(unqday);
    temp_dawn(yrdays,year2do-2002)=dawnhr;
    temp_dusk(yrdays,year2do-2002)=duskhr;
end

dawn_median=floor(nanmedian(temp_dawn,2));
ii=find(dawn_median <8);
dawn_median(ii)=dawn_median(ii-1);

dusk_median=floor(nanmedian(temp_dusk,2));
%fix for dusk data in middle of summer:
ii= dusk_median == 1;
dusk_median(ii)=25; %1 corresponds to first hour of next day :)
ii=find(dusk_median < 20); 
dusk_median(ii)=dusk_median(ii-1); %change value to near neighbor
%%
%save /Volumes/Lab_data/MVCO/FCB/Syn_divrate_model/median_dawn.mat med_dawn med_dusk
save ~/Documents/phyto-division-rate-model/median_dawn.mat dawn_median dusk_median