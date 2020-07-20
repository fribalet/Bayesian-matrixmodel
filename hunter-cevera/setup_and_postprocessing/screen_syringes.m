% little script to walk through record by record?

%first find dates contained within each time/merged file:
year2do=2011;
timepath=['/Volumes/Lab_data/MVCO/FCB/MVCO_Jan2011/data/processed/time/'];
filelist=dir([timepath '*.mat']);

timefile_rec=cell(length(filelist),2);
for q=1:length(filelist)
        
        filename=filelist(q).name;
        eval(['load ' timepath filename])
        disp(['loaded...' filename])
        timefile_rec{q,1}=filename;
        
        %for easier handling, declare time info as temptime:
        if year2do <=2005
            eval(['temptime=' filename(1:end-6) ';'])
        else
            tempvar=whos('FCB*');
            eval(['temptime=' tempvar.name ';'])
        end
        
        timefile_rec{q,2}=[temptime(1,2) temptime(end,3);];
        clear FCB* temptime
end

%% if have dates of interest:

days1=datenum('Jul-10-2011');
days2=datenum('Jul-14-2011');

timestamps=cell2mat(timefile_rec(:,2));
tt=find(timestamps(:,1) < days1 & days2 < timestamps(:,2));
if isempty(tt) %split data chunk?
   tt1=find(timestamps(:,1) < days1);
   tt2=find(days2 < timestamps(:,2));
   tt=[tt1(end);tt2(1)]; %load these two files
end

files2load=timefile_rec{tt,1};
%% mergedwithclass has hourly binned data already...

%display record by record:

for j=150:length(mergedwithclass)
    mwc=mergedwithclass{j};
    
    strec=mwc(1,1);
    endrec=mwc(end,1);
    
    ii=find(FCB1_2011_1time(:,1)==strec-1);
    jj=find(FCB1_2011_1time(:,1)==endrec);
    
    figure(1),clf,loglog(mwc(:,5),mwc(:,2),'.')
    title(['time: ' datestr(FCB1_2011_1time(ii,2)) ' - ' datestr(FCB1_2011_1time(jj,2))])
    figure(2),clf,plot(FCB1_2011_1time(:,2),syrpumpinfo(:,5),'.-'), hold on
    plot(FCB1_2011_1time(ii:jj,2),syrpumpinfo(ii:jj,5),'.-'), datetick
    pause
end