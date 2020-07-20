%Dilution Series Scripts:

%% Individual Bottle Plots:
%T0
ds=2;
eval(['DS_T0C=DS' num2str(ds) '_T0C;'])
eval(['DS_T0B=DS' num2str(ds) '_T0B;'])
eval(['DS_T0A=DS' num2str(ds) '_T0A;'])

eval(['DS_T24C=DS' num2str(ds) '_T24C;'])
eval(['DS_T24B=DS' num2str(ds) '_T24B;'])
eval(['DS_T24A=DS' num2str(ds) '_T24A;'])

%% For the picoeuk data:
ds=2;
eval(['DS_T0C=DS' num2str(ds) '_picoeuk_T0C;'])
eval(['DS_T0B=DS' num2str(ds) '_picoeuk_T0B;'])
eval(['DS_T0A=DS' num2str(ds) '_picoeuk_T0A;'])

eval(['DS_T24C=DS' num2str(ds) '_picoeuk_T24C;'])
eval(['DS_T24B=DS' num2str(ds) '_picoeuk_T24B;'])
eval(['DS_T24A=DS' num2str(ds) '_picoeuk_T24A;'])

%%
eval(['subplot(22' num2str(ds-1) ')'])
h1=plot([0.1 0.3 0.5 0.7 1],(DS_T0C/0.12),'c.','MarkerSize',18);
hold on
plot([0.1 0.3 0.5 0.7 1],(DS_T0B/0.12),'c.','MarkerSize',18)
plot([0.1 0.3 0.5 0.7 1],(DS_T0A/0.12),'c.','MarkerSize',18)

h2=plot([0.1 0.3 0.5 0.7 1],(DS_T24C/0.12),'g.','MarkerSize',18);
plot([0.1 0.3 0.5 0.7 1],(DS_T24B/0.12),'g.','MarkerSize',18)
plot([0.1 0.3 0.5 0.7 1],(DS_T24A/0.12),'g.','MarkerSize',18)

%FCB1 Start and End points:
h3=plot([0.1 1], [Start_End_Points_Diluted_Picoeuks(ds,2) Start_End_Points_Undiluted_Picoeuks(ds,2)],'bo','MarkerSize',10); %beginning points
h4=plot([0.1 1], [Start_End_Points_Diluted_Picoeuks(ds,3) Start_End_Points_Undiluted_Picoeuks(ds,3)],'ro','MarkerSize',10);    %end points

set(gca,'Fontsize',16)
xlabel('Fraction of Unfiltered Seawater','Fontsize',16)
ylabel('Synechococcus Cells per mL','Fontsize',16)
title(['Dilution Series ' num2str(ds)],'Fontsize',16)
% legend([h1 h2 h3 h4],'Time 0H STart Points','Time 24H Endpoints','FCB1 Start Point','FCB1 End Point','Location','NorthWest')

%% growth rates:

%scatter- per bottle growth rate
plot([0.1 0.3 0.5 0.7 1],log(DS_T24C./DS_T0C),'*')
hold on
plot([0.1 0.3 0.5 0.7 1],log(DS_T24B./DS_T0B),'*')
plot([0.1 0.3 0.5 0.7 1],log(DS_T24A./DS_T0A),'*')
%%
figure(3)

hold on
h2=plot([0.1 0.3 0.5 0.7 1],log(nanmean([DS_T24A; DS_T24B; DS_T24C])./nanmean([DS_T0A; DS_T0B; DS_T0C])),'k*--','MarkerSize',16);

%% FCB1 growth rates:

h6=plot([0.1 1],[log(Start_End_Points_Diluted_Picoeuks(ds,3)./Start_End_Points_Diluted_Picoeuks(ds,2)) log(Start_End_Points_Undiluted_Picoeuks(ds,3)./Start_End_Points_Undiluted_Picoeuks(ds,2))],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
%%

set(gca,'Fontsize',16)
xlabel('Fraction of Unfiltered Seawater','Fontsize',16)
ylabel('Growth Rate d^{-1}','Fontsize',16)
title('PicoEuk Diltuion Series Growth Rates','Fontsize',16)
legend([h2 h3 h4 h5 h6],'DS2','DS3','DS4','DS5','FCB1')


%% Fit a line through start points:

figure(2)
plot([0.1 0.3 0.5 0.7 1],(DS_T0A/0.12),'k*','MarkerSize',8);
hold on
plot([0.1 0.3 0.5 0.7 1],(DS_T0B/0.12),'k*','MarkerSize',8)
plot([0.1 0.3 0.5 0.7 1],(DS_T0C/0.12),'k*','MarkerSize',8)

set(gca,'Fontsize',14)
xlabel('Fraction of Unfiltered Seawater','Fontsize',14)
ylabel('Synechococcus Cells per mL','Fontsize',14)
title('Time Zero Points and Best Fit Lines','Fontsize',14)


[R,M,B] = regression([0.1 0.3 0.5 0.7 1],[DS_T0C; DS_T0B; DS_T0A]/.12);
f= @(x) M*x + B;
h1=plot([0 0.1 0.3 0.5 0.7 1],f([0 0.1 0.3 0.5 0.7 1]),'k','Linewidth',1.2);

legend([h4 h3 h2 h1], 'DS2','DS3','DS4','DS5')

%%
h1=plot([0.1 0.3 0.5 0.7 1],(DS3_T0A/0.12),'*','MarkerSize',8);
hold on
plot([0.1 0.3 0.5 0.7 1],(DS3_T0B/0.12),'*','MarkerSize',8)

plot([0.1 0.3 0.5 0.7 1],(DS3_T0C/0.12),'*','MarkerSize',8)


[R,M,B] = regression([0.1 0.3 0.5 0.7 1],[DS3_T0C; DS3_T0B; DS3_T0A]/.12);
f= @(x) M*x + B;
h3=plot([0 0.1 0.3 0.5 0.7 1],f([0 0.1 0.3 0.5 0.7 1]),'--');
legend([h1 h2 h3 h4],'Exp 3 T0', 'Exp 4 T0', 'Best fit line Exp 3' ,'Best fit line Exp 4')