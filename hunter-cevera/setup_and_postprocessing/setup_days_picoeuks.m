
%created 6/27/03 modeled after modelbatch/xiopt4 combination run on LEO-15
%2001 data to save daily size dist files for future model runs, Heidi
%modified 11/4/03 to correct error where same dielstarthr (first one) was
%used for an entire data file; also changed for new time as UTC
%6/24/13-Kristen-changed bead match procedure so that only one beadmatch value is
%used for the entire day. Beadmatches come from smoothed (moving average) mean bead SSC
%using 4 surrounding points to smooth
%Code further modified to use as 'one set' of scripts for all the years
%with user changed inputs for each year. Wrapper script is divmodel_setup.m

%% first, load bead data:
eval(['load ' beadpath 'beadresults.mat']);

%smooth mean of bead SSC with 3-point moving average:

% but first, check for outliers and other bead anomalies:
if plotflag ==1
    figure(6)
    subplot(2,1,1,'replace')
    plot(beadresults(:,1),beadresults(:,13),'.-')
    keyboard
end

[ss is]=sort(beadresults(:,1)); %happens for 2006 where out of time sync?
beadresults=beadresults(is,:);

%known bead outliers:
switch year2do
    case 2003
       ind=find(beadresults(:,13) < 0.9e4); %outlier
       beadresults(ind,13)=NaN;
    case 2009
        ind=find(beadresults(:,13) > 4.05e4); %outlier
        beadresults(ind,13)=NaN; 
    case 2011
        ind=find(beadresults(:,13) > 7e4); %outlier
        beadresults(ind,13)=NaN; 
    case 2013
        ind=find(beadresults(:,13) > 10e4); %outlier
        beadresults(ind,13)=NaN; 
end

sm_bead_avgSSC=mvco_running_average(beadresults(:,1),beadresults(:,13),3,2); %running average smoothing function that takes into account breaks in FCB deployments
%totalrunavg=mvco_running_average(time,data,windowsize,gapsize)
bead_days=floor(beadresults(:,1));

if plotflag ==1
    subplot(2,1,2,'replace')
    hold on
    plot(beadresults(:,1),beadresults(:,13),'.-')
    plot(beadresults(:,1),sm_bead_avgSSC,'.--')
    pause
end

%%
if year2do <= 2005
    filelist = dir([groupedpath '*_*.mat']);
else
    filelist = dir([groupedpath 'FCB*_' num2str(year2do) '*.mat']);
end

%
%Older -fails for some years:
% [junk,ind] = sort(datenum(char(filelist.date)));
% filelist = filelist(ind);

% Get the dates that go along with each file:
%if load straight from filelist, the actual dates will be out of order as
%the way the files line up are in order of FCB#...

groupdates={};
for q=1:length(filelist)
    groupedfile=filelist(q).name;
    eval(['load(''' groupedpath groupedfile ''',''cellresults'')']) %for the date
    ind=find(~isnan(cellresults(:,1)));
    if ~isempty(ind)
        groupdates=[groupdates; {min(cellresults(:,1))} {max(cellresults(:,1))} {groupedfile}];
        %grouplist=[grouplist; repmat(groupedfile,length(ind),1)];
        %filedate=[filedate; cellresults(ind,1)];
    end
    clear cellresults
end

[~, ch_order]=sort(cell2mat(groupdates(:,1))); %sort files into chronological order (chorder is the indexing for this)
filelist=filelist(ch_order);

%
%Load in ligth data:
%eval(['load ' rootpath 'model_code/solar.mat'])
solarpath=fullfile(datapath,'model/');
eval(['load ' solarpath 'solar' num2str(year2do) '.mat']) %QC'd solar for time and nighttime noise

Eall = [date_met Solar];  %date (UTC, matlab), SW radiation
Eall(Eall(:,2)<0,2) = 0;  %no negatives allowed

%volbins = 2.^(-5:0.125:2); %spacing for bins for Syn - log2 base
%volbins = 2.^[-3:0.125:4]; maybe %4.5?
%volbins=2.^[-5:1/6:4.5]; %58 bins - spacing for picoeuks, but too many
%empty small bins

filenum = 1;
filename = filelist(filenum).name;
disp(['loading...' filename])
eval(['load ' groupedpath filename])

%make sure start with non-nan cellresults!
while all(isnan(cellresults(:,1)))
    filenum=filenum+1;
    filename = filelist(filenum).name;
    disp(['loading...' filename])
    eval(['load ' groupedpath filename])
end

if year2do <= 2005
    temp=regexp(filename,'(?<name>[a-z]{2}\d{4}[a-z]{1})(?<num>_\d{1})','names');
    eval(['load ' mergedpath0 temp.name 'merged' temp.num '.mat'])
else
    temp=regexp(filename,'(?<name>FCB\d{1}_\d{4}_\d{1})(?<num>_\d{1})','names');
    eval(['load ' mergedpath0 temp.name 'merged' temp.num '.mat'])
end

%Fixed July 2016 - but before, 2003 was still local time...catch here:
% if year2do == 2003
%     disp('Changing time from local to UTC')
%     if max(cellresults(:,1)) < datenum('Oct-26-2003') %dyalight savings time
%         cellresults(:,1)=cellresults(:,1)+4/24;
%     elseif min(cellresults(:,1)) > datenum('Oct-26-2003') %no daylight savings time
%         cellresults(:,1)=cellresults(:,1)+5/24;
%     elseif  min(cellresults(:,1)) < datenum('Oct-26-2003') & max(cellresults(:,1)) > datenum('Oct-26-2003') %daylight savings time is in this file...
%         ds=find(cellresults(:,1) < datenum('Oct-26-2003'));
%         cellresults(ds,1)=cellresults(:,1)+4/24;
%         cellresults(ds(end)+1:end,1)=cellresults(:,1)+5/24;
%     end
% end


%if ~(strcmp(filename(1:4), 'oc01') | strcmp(filename(1:4), 'oc17') | strcmp(filename(1:4), 'oc19')),
%    cellresults(:,1) = cellresults(:,1) - 5/24;  %convert UTC to local
%end;
%change to all UTC, 11/4/03
%if (strcmp(filename(1:4), 'oc01') | strcmp(filename(1:4), 'oc17') | strcmp(filename(1:4), 'oc19')),
%    cellresults(:,1) = cellresults(:,1) + 4/24;  %convert local to UTC
%end;

daylog_fcb=cell(366,2);


cellclass = 4;%1;  %cell group to model (1 = syn, 4 & 5 = euks)
numfiles = length(filelist);

daylist = floor(min(cellresults(:,1))):floor(max(cellresults(:,1)));
dayind = 1;
day = daylist(dayind);
w=find_yearday(day);

dielstarthr = dawn(dawn(:,1) == day,2);
if isnan(dielstarthr), dielstarthr = 0; end;
if isempty(dielstarthr), dielstarthr = 0; end;
offsetdate = cellresults(:,1) - dielstarthr/24;  %reset so diel start = 0 hour of day

%%
while filenum <= numfiles && day <= floor(max(cellresults(:,1)))  % keep going until past last day in last file + special case of all nan cellresults
    %
    clear cellvol
    mergedpath = mergedpath0;
    dielind = find(floor(offsetdate) == day); %find all hours that are ~24 hours after this dawn
    
    %If have come to the last day in current file, load the next one!
    while day >= floor(max(cellresults(:,1)))-1 && filenum ~= numfiles,  %last day in file, keep partial day and open next file + special case of all nan cellresults
        
        if isempty(dielind) %rare case where second to last day in daylist is empty, so just need to walk back until find data!
            k=0;
            while isempty(dielind) && k <= length(daylist)
                k=k+1;
                dielind = find(floor(offsetdate) == day-k);
            end
            oldcellresults = cellresults(dielind(1):end,:); %keep last day (fix bug from dielind to dielind(1):end 12/11/04)
            oldmergedwithclass = mergedwithclass(dielind(1):end);  %keep last day
            oldbeadmatch = beadmatch(dielind(1):end,:);  %keep last day
        else
            oldcellresults = cellresults(dielind(1):end,:); %keep last day (fix bug from dielind to dielind(1):end 12/11/04)
            oldmergedwithclass = mergedwithclass(dielind(1):end);  %keep last day
            oldbeadmatch = beadmatch(dielind(1):end,:);  %keep last day
        end
        
        filenum = filenum + 1;
        filename = filelist(filenum).name;  %next file
        disp(['loading...' filename])
        eval(['load ' groupedpath filename])
        
        if year2do <= 2005
            temp=regexp(filename,'(?<name>[a-z]{2}\d{4}[a-z]{1})(?<num>_\d{1})','names');
            eval(['load ' mergedpath0 temp.name 'merged' temp.num '.mat'])
        else
            temp=regexp(filename,'(?<name>FCB\d{1}_\d{4}_\d{1})(?<num>_\d{1})','names');
            eval(['load ' mergedpath0 temp.name 'merged' temp.num '.mat'])
        end
        
        %older fix, no longer necessary since July 2016
%         if year2do == 2003
%             disp('Changing time from local to UTC')
%             if max(cellresults(:,1)) < datenum('Oct-26-2003') %dyalight savings time
%                 cellresults(:,1)=cellresults(:,1)+4/24;
%             elseif min(cellresults(:,1)) > datenum('Oct-26-2003') %no daylight savings time
%                 cellresults(:,1)=cellresults(:,1)+5/24;
%             elseif  min(cellresults(:,1)) < datenum('Oct-26-2003') & max(cellresults(:,1)) > datenum('Oct-26-2003') %daylight savings time is in this file...
%                 ds=find(cellresults(:,1) < datenum('Oct-26-2003'));
%                 cellresults(ds,1)=cellresults(:,1)+4/24;
%                 cellresults(ds(end)+1:end,1)=cellresults(:,1)+5/24;
%             end
%         end
        
        cellresults = [oldcellresults; cellresults];  %add last day to start of next file
        mergedwithclass = [oldmergedwithclass mergedwithclass];
        
        if isnan(dielstarthr), dielstarthr = 0; end;
        if isempty(dielstarthr), dielstarthr = 0; end;
        
        offsetdate = cellresults(:,1) - dielstarthr/24;  %recalculate all dates and indices
        dielind = find(floor(offsetdate) == day);
        dayind = 1; %restart at beginning of daylist
        clear old*
    end;  %if day = floor(cellresults(end,1)),
   %
    disp(['day: ' num2str(day)])
    daylist = floor(min(cellresults(:,1))):floor(max(cellresults(:,1)));
    
    if ~isempty(dielind) && dielind(end) < size(offsetdate,1),
        dielind = [dielind; dielind(end)+1];
        dielhr = floor((offsetdate(dielind) - day)*24);
        dielind = dielind(dielhr <= 24);
        dielhr = floor((offsetdate(dielind) - day)*24);
    else   %case where no data points for day or last date
        dielhr = NaN;
        daylog_fcb(w,:)=[{day} {'no fcb data'}];
    end;
    %check to see if have beginning hour 0:
    if dielhr(1) == 1,  %repeat first hr 1 for hr 0
        dielind = [dielind(1); dielind];
        dielhr = [0; dielhr];
        disp('missing first hour')
        daylog_fcb(w,:)=[{day} {'missing first hour'}];
    end;
    %check to see if have ending hour 24:
    if dielhr(end) == 23,  %repeat hr 23 as hr 24
        dielind = [dielind; dielind(end)];
        dielhr = [dielhr; 24];
        disp('missing last hour')
        if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day} {'missing last hour'}];
        else daylog_fcb{w,2}=[daylog_fcb{w,2} '; missing first hour']; end
    end;
    %
    %Is a full day or last 19 hours present?
    if dielhr(1) <= 5 && dielhr(end) == 24 && isempty(find(diff(dielhr) > 2, 1)),  %full day ==> proceed, otherwise just go to next day
        %if dielhr(1) == 0 & dielhr(end) == 24 & isempty(find(diff(dielhr) > 2)),  %full day ==> proceed, otherwise just go to next day
        
        cellresultsfordiel = cellresults(dielind,:);
        
        for count = 1:length(dielind),  %get all cell data
            temp = mergedwithclass{dielind(count)};
            ind = find(temp(:,end) == cellclass);
            
            %find matching bead data for normalization:
            beadind=find(bead_days==day);
            if ~isempty(beadind) %found beads
                if length(beadind)==1
                    beadvol=sm_bead_avgSSC(beadind); %use smoothed bead data
                else
                    beadvol=mean(sm_bead_avgSSC(beadind)); %if two events were measured, avg over them
                end
                
            else %missing day in bead data- use bead values close to that day or average around that day:
                bi=find(bead_days==day-1);
                bii=find(bead_days==day+1);
                if isempty(bi) && isempty(bii)
                    bi=find(bead_days==day-2);
                    bii=find(bead_days==day+2);
                end
                if ~isempty(bi) && ~isempty(bii)
                    beadvol=mean([sm_bead_avgSSC(bi(end)) sm_bead_avgSSC(bii(1))]);
                elseif isempty(bi) & ~isempty(bii)
                    beadvol=mean(sm_bead_avgSSC(bii(1)));
                elseif isempty(bii) & ~isempty(bi)
                    beadvol=mean(sm_bead_avgSSC(bi(end)));
                end
            end
            
            cellvol{count} = cytosub_SSC2vol(temp(ind,5)./beadvol);  %SSC bu values converted to volume
            %cellvol{count} = cytosub_SSC2vol(temp(ind,5)./beadmatch(dielind(count),5));  %SSC bu values converted to volume
        end; %for count
        %
        temp = find(diff(dielhr) == 0);  %case where more than one sample in same hr (only at quick acq restart)
        while ~isempty(temp), %every other one matches
            temp = temp(1); %deal with first one first...loop to get subsequent
            cellvol{temp} = [cellvol{temp}; cellvol{temp+1}]; %add together cells measured in same hour
            cellvol = cellvol([1:temp, temp+2:end]); %omit second in same hour
            %                cellresultsfordiel(temp,:) = mean(cellresultsfordiel(temp:temp+1,:));
            cellresultsfordiel(temp,:) = [mean(cellresultsfordiel(temp:temp+1,1)) cellresultsfordiel(temp,2)+cellresultsfordiel(temp+1,2) mean(cellresultsfordiel(temp:temp+1,3))];  %add together acq time...
            cellresultsfordiel = cellresultsfordiel([1:temp, temp+2:end],:);
            dielhr = dielhr([1:temp, temp+2:end]);
            dielind = dielind([1:temp, temp+2:end]);
            temp = find(diff(dielhr) == 0);  %check again
            disp('more than one sample in same hr...')
            if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day}  {'more than one sample in hour'}];
            else daylog_fcb{w,2}=[daylog_fcb{w,2} '; more than one sample in hour']; end
        end;
        
        temp = find(diff(dielhr) == 2);  %case where one hour missing
        while ~isempty(temp),  %average across 1 skipped hour
            temp = temp(1);
            tempvol{1} = [cellvol{temp}; cellvol{temp+1}]; %add together cells measured in surrounding hours
            cellvol = [cellvol([1:temp]) tempvol cellvol([temp+1:end])];
            tempcellres = mean(cellresultsfordiel([temp, temp+1],:));
            tempcellres(2) = cellresultsfordiel(temp,2) + cellresultsfordiel(temp+1,2);  %add together acq times
            cellresultsfordiel = [cellresultsfordiel(1:temp,:); tempcellres; cellresultsfordiel(temp+1:end,:)];
            dielhr = [dielhr(1:temp); dielhr(temp)+1; dielhr(temp+1:end)];
            dielind = [dielind(1:temp); dielind(temp)+1; dielind(temp+1:end)];
            temp = find(diff(dielhr) == 2); %look for another...
            disp('  averaged over one missing hour')
            if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day}  {'avg over missing hour'}];
            else daylog_fcb{w,2}=[daylog_fcb{w,2} '; avg over missing hour']; end
        end;
        clear temp count
        
        %match corresponding light data:
        Eind = find(floor(Eall(:,1) - dielstarthr/24) == day);  %offset Edata to dawn
        if ~isempty(Eind),
            %for rare case when end of light data occurs during a day:
            if Eind(end)+1 > length(Eall)
                Eind=Eind(1):length(Eall);
            else
                Eind = [Eind; Eind(end)+1];
            end            
            Edata = Eall(Eind,:);
            if (Eall(Eind(end),1) - Eall(Eind(end-1))) > 2/24, %make sure end of day does not interpolate to middle of next light
                Edata(end,2) = 0;
            end;
            Edata(:,1) = (Edata(:,1) - day)*24 - dielstarthr;  %Radiation data for hours from 0 to 24
        else
            Edata = [NaN NaN];
        end;
        clear Eind
        
        Vhists = zeros(length(volbins),length(cellvol));
        N_dist = Vhists;
        
        for i = 1:length(cellvol) %calculate observed size distributions
            N_dist(:,i) = histc(cellvol{i}, volbins);  %size distribution
            Vhists(:,i) = N_dist(:,i)/sum(N_dist(:,i));    %normalized size distribution
        end
        %cellsperml = [sum(N_dist)./cellresultsfordiel(:,2)'./cellresultsfordiel(:,3)'];
        cellsperml = [sum(N_dist)./cellresultsfordiel(:,3)'];
        
        if length(cellvol)==25
            if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day}  {'full day!'}];
            else daylog_fcb{w,2}=[daylog_fcb{w,2} '; full day!']; end
        elseif length(cellvol) <25
            if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day}  {'partial day!'}];
            else daylog_fcb{w,2}=[daylog_fcb{w,2} '; partial day!']; end
        end
        
        if dielstarthr == 0,
            disp(['...skip; Missing Edata '])
            if isempty(daylog_fcb{w,2}), daylog_fcb(w,:)=[{day}  {'missing Edata'}];
            else daylog_fcb{w,2}=[daylog_fcb{w,2} '; missing Edata']; end
        else
 %           eval(['save ' modelpath 'day' num2str(day) 'data volbins Edata Vhists N_dist cellsperml dielstarthr'])  %save data for rerunning batch later
        end;
    else
        disp(['...skip; dielhr: ' num2str(dielhr')])
    end; %if dielhr(1) <= 5 & dielhr(end) == 25,  %full day ==> proceed
    %
    if filenum == numfiles && day >= daylist(end), %case of last file and last day (untested as of 2/16/02)
        day = day + 1;
        w=find_yearday(day);
        %        dielstarthr = dawn(find(dawn(:,1) == day),2);
        %        offsetdate = cellresults(:,1) - dielstarthr/24;  %reset so diel start = 0 hour of day
    else  %normal case
        dayind = dayind + 1;  %go to next day
        day = daylist(dayind);
        w=find_yearday(day);
        dielstarthr = dawn(find(dawn(:,1) == day),2);
        if isnan(dielstarthr), dielstarthr = 0; end;
        if isempty(dielstarthr), dielstarthr = 0; end;
        %inserted here 11/4/03 to fix previous error (dielstarthr = first value for entire file...)
        offsetdate = cellresults(:,1) - dielstarthr/24;  %reset so diel start = 0 hour of day
    end;
    
end;     %while filenum <= numfiles & day <= floor(cellresults(end,1)),  % keep going until past last day in last file

%Save log for reference:
emptydays=find(cellfun('isempty',daylog_fcb(:,1))==1);
daylog_fcb(emptydays,1)=num2cell(datenum(repmat(['1-0-' num2str(year2do)],size(emptydays,1),1))+emptydays); %date
daylog_fcb(emptydays,2)=cellstr(repmat('no fcb data',size(emptydays,1),1));
eval(['save ' modelpath 'setupdays_log.mat daylog_fcb'])