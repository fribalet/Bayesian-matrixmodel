%for year that had a few days to rerun:
clear
year=2014; yearlabel='Jan';

eval(['cd /Volumes/Lab_data/MVCO/FCB/MVCO_' yearlabel num2str(year) '/model/output_July2016/'])
eval(['load(''mvco_14par_dmn_' num2str(year) 'redoC.mat'')'])
modelresultsR=modelresults; allmodelrunsR=allmodelruns;
eval(['load(''mvco_14par_dmn_' num2str(year) '.mat'')'])
jj=find(modelresultsR(:,1)~=0)
[modelresultsR(jj,16)-modelresults(jj,16) modelresultsR(jj,17) modelresults(jj,17)]
%% take a look at the modelfits:

for q=1:length(jj)
    tempR=allmodelrunsR{jj(q),1}; temp=allmodelruns{jj(q),1};
    [ss is]=sort(temp(:,15)); [ss ir]=sort(tempR(:,15));
    subplot(1,2,1,'replace'), plot(1:length(tempR),tempR(ir,16),'.-'), hold on, plot(1:length(temp),temp(is,16),'.-')
    subplot(1,2,2,'replace'), plot(1:length(tempR),tempR(ir,15),'.-'), hold on, plot(1:length(temp),temp(is,15),'.-')
    title(['Day: ' datestr(modelresults(jj(q),1)) ' ' num2str(modelresults(jj(q),1))])
    pause
end

%% then if still want to replace:
eval(['system(''mv mvco_14par_dmn_' num2str(year) '.mat mvco_14par_dmn_' num2str(year) '_original.mat'')']) %rename the original runs so won't overwirte
%%

for i=1:length(jj)
    if modelresultsR(jj(i),16) < modelresults(jj(i),16)
        modelresults(jj(i),:)=modelresultsR(jj(i),:);
    end
    %combine all model runs: 
    allmodelruns{jj(i),1}=[allmodelruns{jj(i),1}; allmodelrunsR{jj(i),1}];
    allmodelruns{jj(i),2}=[allmodelruns{jj(i),2}; allmodelrunsR{jj(i),2}];
end

eval(['save mvco_14par_dmn_' num2str(year) '.mat modelresults allmodelruns'])