%Dilution series summary plots:
ds=11;
    if ds<=5 || ds==8 %1 day experiments
        eval(['load DS' num2str(ds) '_rawdata'])
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;']) 
        ds_exp_plot(ds,dilutions,DS_T0,DS_picoeuks_T0,DS_T24,DS_picoeuks_T24)
    else %2-day experiments
        eval(['load DS' num2str(ds) '_rawdata'])
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_T48= DS' num2str(ds) '_T48;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;'])
        eval(['DS_picoeuks_T48= DS' num2str(ds) '_picoeuks_T48;'])
        ds_exp_plot(ds,dilutions,DS_T0,DS_picoeuks_T0,DS_T24,DS_picoeuks_T24,DS_T48,DS_picoeuks_T48)
    end
%%
for ds=2:11
    eval(['load DS' num2str(ds) '_rawdata'])
end
%% Across all DS Exp plots:
color=hsv(10);
c=0;

for ds=2:11
    c=c+1;
    if ds<=5 || ds==8 %1 day experiments
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;'])       
    else %2-day experiments
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_T48= DS' num2str(ds) '_T48;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;'])
        eval(['DS_picoeuks_T48= DS' num2str(ds) '_picoeuks_T48;'])
    end
    [n m]=size(DS_T0);
    
    if ds <=6
        subplot(121)
        hold on
        plot(nanmean(DS_T0),ds*ones(1,m),'*:','Color',color(c,:))
        plot(nanmean(DS_picoeuks_T0),ds*ones(1,m)+0.2,'o:','Color',color(c,:))
        if ds==6
            plot(nanmean(DS_T24),ds*ones(1,m)+0.3,'*:','Color',color(c,:))
            plot(nanmean(DS_picoeuks_T24),ds*ones(1,m)+0.4,'o:','Color',color(c,:))
        end
    else
        subplot(122)
        hold on
        plot(nanmean(DS_T0),ds*ones(1,m),'*:','Color',color(c,:))
        plot(nanmean(DS_picoeuks_T0),ds*ones(1,m)+0.2,'o:','Color',color(c,:))
        if ds~=8
        plot(nanmean(DS_picoeuks_T24),ds*ones(1,m)+0.4,'o:','Color',color(c,:))
        plot(nanmean(DS_T24),ds*ones(1,m)+0.3,'*:','Color',color(c,:))
        end
    end
end
subplot(121)
title('Summer Dilution Experiments','Fontsize',14)
xlabel('Cell Concentration (per mL)','Fontsize',14)
ylabel('Dilution Series Exp #','Fontsize',14)
set(gca,'Fontsize',14)

subplot(122)
title('Fall Dilution Experiments','Fontsize',14)
xlabel('Cell Concentration (per mL)','Fontsize',14)
ylabel('Dilution Series Exp #','Fontsize',14)
set(gca,'Fontsize',14)
%% By Syn number:
color=hsv(10);
c=0;
H1=[];
H2=[];
titles1={};
titles2={};

for ds=2:11
    c=c+1;
    if ds<=5 || ds==8 %1 day experiments
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;']) 
        
        netdiff1=(log(nanmean((DS_picoeuks_T24))./nanmean((DS_picoeuks_T0)))); netdiff1=netdiff1(1)-netdiff1(end);
    else %2-day experiments
        eval(['DS_T0= DS' num2str(ds) '_T0;'])
        eval(['DS_T24= DS' num2str(ds) '_T24;'])
        eval(['DS_T48= DS' num2str(ds) '_T48;'])
        eval(['DS_picoeuks_T0= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_picoeuks_T24= DS' num2str(ds) '_picoeuks_T24;'])
        eval(['DS_picoeuks_T48= DS' num2str(ds) '_picoeuks_T48;'])
        
    netdiff1=(log(nanmean((DS_picoeuks_T24))./nanmean((DS_picoeuks_T0)))); netdiff1=netdiff1(1)-netdiff1(end);
    netdiff2=(log(nanmean((DS_picoeuks_T48))./nanmean((DS_picoeuks_T24)))); netdiff2=netdiff2(1)-netdiff2(end);
    end
    
    if ds <=6
%         subplot(121) %figure(1)
        hold on
%         h1=plot(min(nanmean(DS_T24)),netdiff1,'*:','Color',color(c,:));
        h1=plot(min(nanmean(DS_picoeuks_T24)),netdiff1,'*:','Color',color(c,:));
        H1=[H1; h1];
        titles1=[titles1; num2str(ds)];
        if ds==6
             plot(min(nanmean(DS_picoeuks_T48)),netdiff2,'s:','Color',color(c,:))
        end
    else
%         subplot(122) %figure(2)
        hold on
         h2=plot(min(nanmean(DS_picoeuks_T24)),netdiff1,'*:','Color',color(c,:));
        H2=[H2; h2];
        titles2=[titles2; num2str(ds)];
        if ds~=8
         plot(min(nanmean(DS_picoeuks_T48)),netdiff2,'s:','Color',color(c,:))
        end
    end
end
%    
% figure(1)
% legend(H1,titles1)
% figure(2)
% legend(H2,titles2)
legend([H1; H2], [titles1; titles2])
%% Grazing rate plots
clf
record=[];
for q=1:length(ds_exp_mu)
    ds=ds_exp_mu(q,2);
    if ds_exp_mu(q,3)==1
        eval(['load DS' num2str(ds) '_rawdata'])
        eval(['DS_T1= DS' num2str(ds) '_T0;'])
        eval(['DS_T2= DS' num2str(ds) '_T24;'])
    else 
        eval(['DS_T1= DS' num2str(ds) '_T24;'])
        eval(['DS_T2= DS' num2str(ds) '_T48;'])
    end
    if ~isnan(ds_exp_mu(q,4))
        grazing_rate=ds_exp_mu(q,4)-log(nanmean((DS_T2(:,end)))./nanmean((DS_T1(:,end))));
    else
        grazing_rate=ds_exp_mu(q,5)-log(nanmean((DS_T2(:,end)))./nanmean((DS_T1(:,end))));
    end
    record=[record; ds nanmean(DS_T1(:,end-1)) grazing_rate];
end

%%
color=hsv(4);
i=0;
for j=4:2:10
    i=i+1;
    plot(record(j,2),record(j,3),'s','Color',color(i,:))
    hold on
end