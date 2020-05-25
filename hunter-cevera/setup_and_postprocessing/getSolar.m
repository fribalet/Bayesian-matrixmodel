if year2do == 2016
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99B.txt'',''ascii'');'])
elseif year2do == 2017
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_sB.C99.txt'',''ascii'');'])
else
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99'',''ascii'');'])
end

yd_met=X(:,2);
Hour = X(:,5);
Solar=X(:,14);  %Solar_campmt_median
% Year = X(:,1);
% Day = X(:,4);

clear X

%fix for duplicated data, sometimes interspersed with good data?
[tt, it]=sort(yd_met);
yd_met=yd_met(it);
Solar=Solar(it);
Hour=Hour(it);

[tt2, iu]=unique(yd_met);
yd_met=yd_met(iu);
Solar=Solar(iu);
Hour=Hour(iu);

t = find(diff(yd_met) < 0);
fprintf('Found "out of order" time for %1.0f instances!\n',length(t))
yd_met(t+1) = NaN;
Solar(t+1) = NaN;
Hour(t+1) = NaN;
% Day(t+1) = NaN;
% Year(t+1) = NaN;

%exclude nan's:
ind = find(~isnan(Solar) & ~isnan(yd_met));
yd_met = yd_met(ind);
Solar = Solar(ind);
Hour = Hour(ind);

date_met = yd_met + datenum(['1-0-' num2str(year2do)]);  %should be UTC
dawnlevel = 5;   %threshold for light

Solar0=Solar; %copies....
date_met0=date_met;
%% Load in good estimates of where dawn/dusk is to detect gaps:

%constructed from raw_data with get_average_dawn.m, 
%load expected UTC values of dawn and dusk over 13 years :)
if ~isempty(strfind(computer,'WIN'))
    load \\sosiknas1\lab_data\mvco\FCB\Syn_divrate_model\median_dawn.mat
else %temporary locations:
    load /Users/kristenhunter-cevera/phyto-division-rate-model/setup_and_postprocessing/median_dawn.mat  
end 


%% find dawn and dusk of each day and detect gaps:

days=(datenum(['1-1-' num2str(year2do)]):datenum(['12-31-' num2str(year2do)]))';
duskhr=nan(length(days),1);
dawnhr=nan(length(days),1);
lighttime_gaps=nan(length(days),1);

for j=1:length(days)
    
    ind = find(floor(date_met - 5/24) == days(j));
    if ~isempty(ind)
        
        ind2=find(date_met(ind) >= days(j)+(dawn_median(j)-1)/24 & date_met(ind) <= days(j)+(dusk_median(j)+1)/24); %examine light only within the expected dawn-dusk range (sometimes renegard nightime noise!)
        ind3 = find(Solar(ind(ind2)) > dawnlevel); %is there light data?
        
        if ~isempty(ind3) %good, we have light data!
            
            tempdawn=Hour(ind(ind2(ind3(1)))) - 1; %1 h before dawnlevel?
            tempdusk=Hour(ind(ind2(ind3(end)))) + 1; %next hour after end of light
            if tempdusk < 20, tempdusk=tempdusk+24; end %i.e. for hours 1,2 - this is the next day...
            
            %check for gaps around dawn and dusk:           
            if abs(dawn_median(j)-tempdawn) <= 1
                dawnhr(j)=tempdawn; %good, accept this dawn!            
            elseif abs(dawn_median(j)-tempdawn) <= 3 %close...small gap, pad with a datapoint for later interpolation!
                 date_met = [date_met; days(j)+(dawn_median(j)-.001)/24];
                 Solar = [Solar; 0];
                 dawnhr(j)=dawn_median(j);
            else
                lighttime_gaps(j)=0; %gap around dawn
            end
                        
            if abs(dusk_median(j)-tempdusk) <= 1
                duskhr(j)=tempdusk; %good, accept this dusk!            
            elseif abs(dusk_median(j)-tempdusk) <= 3 %close...small gap, pad with a datapoint for later interpolation!
                 date_met = [date_met; days(j)+(dusk_median(j)+.001)/24];
                 Solar = [Solar; 0];
                 duskhr(j)=dusk_median(j);
            else
                lighttime_gaps(j)=2; %gap around dusk
            end
            
            %now check for any gaps during light time hours:
            if ~isnan(dawnhr(j)) && ~isnan(duskhr(j))
                if any(find(diff(date_met(ind(ind2))) >= 2/24));
                    lighttime_gaps(j) = 1;
                end;
            end
            
            % a quick double check on days 
            %QC plotting if needed:
%             if ~isnan(lighttime_gaps(j))
%                 hold on
%                 plot(date_met0(ind),Solar0(ind),'ko')
%                 plot(date_met0(ind(ind2)),Solar0(ind(ind2)),'r*')              
%                 keyboard
%             end
            
            
        end
    end
end

%in case added any points:
[~,ii]=sort(date_met);
date_met=date_met(ii);
Solar=Solar(ii);


%% %%%%%%% identify missing days, and for years 2005-2007, 2010-2013, check if buoy data can be substituted! %%%%%%%%%%%%%%%%%%%%%%%%

% For most of 2010, need to use nantucket buoy data!
%Buoy data is formatted as:
%#YY  MM DD hh mm  SRAD1  SWRAD  LWRAD
%#yr  mo dy hr mn   w/m2   w/m2   w/m2
%but appears to use different instruments for different years...

if ismember(year2do,[2005:2007 2010:2013]) && buoy_flag==1; %meaning, yes - you'd like to splice in buoy data
    
    disp('Checking Nantucket Buoy data for any missing gaps in MVCO record...')
    filename = fullfile(envpath,'Nantucket_buoy44008_lightdata',['44008r' num2str(year2do)]);
    startRow = 3;
    formatSpec = '%4f%3f%3f%3f%3f%7f%7f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    buoy_data=cell2mat(dataArray(:,1:8));
    buoy_date=datenum([buoy_data(:,1:5) zeros(length(buoy_data(:,1)),1)]);
    buoy_date=buoy_date-((2/3)/24); %seems to better align with MVCO data if shifted
    %data is reported for ending hour, but not sure when in the hour data was taken!
    clearvars filename startRow formatSpec fileID dataArray ans;  
    
    if ismember(year2do,2010:2013)
        sw_rad=buoy_data(:,7);       
    elseif ismember(year2do,2005:2007)
        sw_rad=buoy_data(:,6);
    end
          
    % splice in nantucket buoy data:
    % go through each day and see when data is available - then splice together:
    missing_days=find(~isnan(lighttime_gaps));
    buoy_added=[];
    buoy_days=nan(length(days),1);
    mvco_to_remove=[];
    
    for q=1:length(missing_days) %check if buoy data has them!
        tempday=missing_days(q)+datenum(['1-0-' num2str(year2do)]);     
        qq=find(buoy_date-5/24 >=tempday & buoy_date-5/24 <tempday+1); %local frame of reference!
        
        if ~isempty(qq) && length(qq) >=22 %buoy data is available and seems complete!            
            %add data from buoy            
            buoy_added=[buoy_added; buoy_date(qq) sw_rad(qq)];   
            
            %add to dawn:
            ind2 = find(sw_rad(qq) > dawnlevel);
            
            tempdawnhr=buoy_data(qq(ind2(1)),4)-1-1; %1 h before dawnlevel? (-1 is for offset of data recording)
            if abs(tempdawnhr-dawn_median) > 1, tempdawnhr = dawn_median(missing_days(q)); end
            tempduskhr=buoy_data(qq(ind2(end)),4)+1-1; %next hour after end of light plus offset
            if abs(tempduskhr-dusk_median) > 1, tempduskhr = dusk_median(missing_days(q)); end
            
            dawnhr(missing_days(q)) = tempdawnhr;  
            duskhr(missing_days(q)) = tempduskhr;  
            buoy_days(missing_days(q)) = 1;
            
            %data to remove from MVCO track:
            mm=find(date_met >= buoy_date(qq(1)) & date_met <= buoy_date(qq(end)));
            mvco_to_remove=[mvco_to_remove; mm];
        end   
    end
   
    % now add in and remove data - yikes!  
    [ind]=setdiff(1:length(date_met),mvco_to_remove); %mvco data to keep, that won't be replaced
    date_met_temp=[date_met(ind); buoy_added(:,1)];
    solar_temp=[Solar(ind); buoy_added(:,2)];  
    [~, is]= sort(date_met_temp);
    date_met=date_met_temp(is);
    Solar=solar_temp(is);
    
end

dawn=[days dawnhr];
ii=find(isnan(buoy_days));
jj=find(~isnan(lighttime_gaps(ii)));
dawn(ii(jj),2)=NaN; %not a good day

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  check for nighttime noise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

date_met_local=date_met-5/24; %for easier handling

for j=1:length(days)
    
    day=days(j);
    ii=find(date_met_local >= day & date_met_local <= day+1);
    
    %find any light levels during expected nightime hours:
    nn=find(date_met_local(ii) <= day+dawn_median(j)/24-5/24 | date_met_local(ii) >= day+dusk_median(j)/24-5/24); %min normal dawn and max normal dusk

    if ~isempty(nn)
        if any(Solar(ii(nn)) > 20)
            disp(['Found "light" readings during dark period for day: ' num2str(day) ': ' datestr(day) ' correcting...'])
            Solar(ii(nn)) = 0;            
            %             plot(date_met_local(ii(nn)),Solar(ii(nn)),'.-')
            %             xlim([day-3 day+3])
            %             datetick('x','mm dd','keeplimits')
            %             pause
            
        end
    end
end

%negative noise, simply remove:
jj=find(Solar < 0);
disp(['replacing ' num2str(length(jj)) ' negative values with 0'])
Solar(jj)=0;

%% %%%%%%%%%%%%%%%%%%%%%%   AND PLOT!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if solarplotflag
    clf
    %colorblock expected nightime:
    for j=1:366
        day=datenum(['1-0-' num2str(year2do)])+j;
        %color in night:
        f1=fill([day-1+dusk_median(j)/24; day+dawn_median(j)/24; day+dawn_median(j)/24; day-1+dusk_median(j)/24],[0 0 1000 1000],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    
    %plot the data!
    h1=plot(date_met,Solar,'-','color',[0.0265 0.6137 0.8135]);
    set(gca, 'layer', 'top')
    h2=plot(date_met0,Solar0,'.','color',[0 0.4 0.8]);
    
    if exist('buoy_added','var')
        hold on
        h3=plot(buoy_added(:,1),buoy_added(:,2),'o','linewidth',2,'markersize',4,'color',[0 0 0.7]);
        legend([h1(1); h2(1); h3(1)],'Final data','MVCO data','Original buoy data','location','NorthOutside')
        title([num2str(year2do) ' Solar data: ' num2str(length(find(~isnan(buoy_days)))) ' days added from buoy data'])
    else
        legend([h1(1); h2(1)],'Final data','MVCO data','location','NorthOutside')
        title(num2str(year2do))
    end
    
    %add dawn lines:
    for i=1:length(dawn)
        if isnan(dawn(i,2))
           plot(dawn(i,1)+(dawnhr(i)+12)/24,800,'x','linewidth',2,'color',[0 0 0]) %first matlab default color: [0 0.5 0.8]
        else
           line([dawn(i,1)+dawn(i,2)/24 dawn(i,1)+dawn(i,2)/24],[0 1000],'linewidth',2,'color',[0.8 0.5 0]) %first matlab default color: [0 0.5 0.8]
        end
    end
    
    
end

xlim([datenum(['1-1-' num2str(year2do)]) datenum(['12-31-' num2str(year2do)])])
ylim([-10 max(Solar)])
datetick('x','keeplimits')




%% %SAVE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('buoy_added','var')
    eval(['save ' solarsavepath 'solar' num2str(year2do) '_w_buoy.mat Solar date_met dawn buoy_added'])
else
    eval(['save ' solarsavepath 'solar' num2str(year2do) '.mat Solar date_met dawn'])
end



%% screen days, one by one:
%  for count = 1:length(unqday),
%     ind = find(floor(date_met) == unqday(count));
%     figure(1), clf
%     plot(yd_met(ind), Solar(ind), '.-')
%     hold on
%     line([dawnhr(count) dawnhr(count)]/24+unqday(count)-datenum(['1-0-' num2str(year2do)]), [0 1000], 'color', 'r')
%     line([duskhr(count) duskhr(count)]/24+unqday(count)-datenum(['1-0-' num2str(year2do)]), [0 1000], 'color', 'g')
%     title([num2str(unqday(count)) ' ; ' datestr(unqday(count))])
%     pause
%  end;
%