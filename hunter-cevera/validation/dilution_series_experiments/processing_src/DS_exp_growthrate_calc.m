%Dilution Experiment Linear and 3-point growth fits:

% DS_T2=DS9_T48*(1/0.36);
% DS_T1=DS9_T24*(1/0.36);

DS_T2=DS11_picoeuks_T48;
DS_T1=DS11_picoeuks_T24;

% DS_T1=[DS_T1; nan(1,5)]; %only for DS3
% DS_T1=DS_T1(1:3,:); %only for DS5

% DS_T2=DS7_T48;
% DS_T1=DS7_T24;

%  dilution=  [0.1 0.3 0.5 0.7 1];
 dilution=[0.2 0.4 0.6 1];
% dilution=[0.1 0.3 0.6 1];
% dilution=dilutions

%straight linear fit:
Y=log((DS_T2)./(DS_T1));
% exclude an outlier point for DS10 day1:
% Y(3,1)=NaN;

x=repmat(dilution,3,1);
x=x(:);
X=[ones(length(x),1) x]; %add coefficent to do regression

[blin,~,~,~,stats] = regress(Y(:),X);
clear bint rint r
%b(1) is the intercept, b(2) is the slope!
%stats: r2=stats(1), pvalue=stats(2)
clf
plot(dilution,log((DS_T2)./(DS_T1)),'k*')
hold on
plot([0 dilution],blin(1)+blin(2)*[0 dilution],'b--')

%%
%linear fit with last two points (aka 3 point method)
Y=log((DS_T2(:,1:2))./(DS_T1(:,1:2)));
%exclude an outlier point for DS10 day1:
% Y(3,1)=NaN;
x=repmat(dilution(1:2),3,1);
x=x(:);
X=[ones(length(x),1) x]; %add coefficent to do regression

[b3,~,~,~,stats] = regress(Y(:),X);
clear bint rint r
%b(1) is the intercept, b(2) is the slope!
%stats: r2=stats(1), pvalue=stats(2)
% plot(dilution,log((DS_T2)./(DS_T1)),'k*')
% hold on
plot([0 dilution],b3(1)+b3(2)*[0 dilution],'g')
ylim([min(min(log((DS_T2)./(DS_T1))))-0.2 b3(1)+0.2])
line([0 1], [0 0],'LineStyle','--','Color','k')
set(gca,'Fontsize',16)
xlabel('Fraction Unfiltered Seawater','Fontsize',14)
ylabel('Net Growth Rate','Fontsize',14)
%%
day=
ds_exp_mu=[ds_exp_mu; 11 2 b3(1) blin(1)]
%%

ds_exp_mu_picoeuks=[ds_exp_mu_picoeuks; 11 2 b3(1) blin(1)]
%%
figure

IND=[];
subplot(121)
hold on
for j=1:length(ds_exp_mu)
    ind=find(ds_exp_mu(j)==MR_12params_tank3(:,1));
    IND=[IND; ind];
    h1=plot(ds_exp_mu(j,4),MR_12params_tank3(ind,15),'g*');
    h2=plot(ds_exp_mu(j,4),MR_12params_tank6(ind,15),'gs');
end

line([0 1],[0 1])
set(gca,'Fontsize',14)
legend([h1(1) h2(1)], '2-point fit Diluted','2-point fit Undiluted')
xlabel('Dilution Series Exp Growth Rate','Fontsize',14)
ylabel('Model Growth Rate','Fontsize',14)

subplot(122)
h3=plot(ds_exp_mu(:,5),MR_12params_tank3(IND,15),'b*');
hold on
h4=plot(ds_exp_mu(:,5),MR_12params_tank6(IND,15),'bs');

legend([h3(1) h4(1)],'Linear Regression fit Diluted','Linear Regression Undiluted')
xlabel('Dilution Series Exp Growth Rate','Fontsize',14)
line([0 1],[0 1])
ylabel('Model Growth Rate','Fontsize',14)

%%
figure
subplot(121)
plot(MR_12params_tank3(IND,15), MR_12params_tank6(IND,15),'ko')
line([0 1],[0 1])
set(gca,'Fontsize',14)
xlabel('Diluted Growth Rate','Fontsize',14)
ylabel('Undiluted Growth Rate','Fontsize',14)

%% Either 3-point or linear fit:

subplot(1,2,2,'replace')
hold on

for j=1:length(ds_exp_mu)
    if ~isnan(ds_exp_mu(j,4))
    h1=plot(ds_exp_mu(j,4),MR_12params_tank3(IND(j),15),'g*');
    h2=plot(ds_exp_mu(j,4),MR_12params_tank6(IND(j),15),'gs');
    else
    h3=plot(ds_exp_mu(j,5),MR_12params_tank3(IND(j),15),'*');
    h4=plot(ds_exp_mu(j,5),MR_12params_tank6(IND(j),15),'s');
    end
end
%%
line([0 1],[0 1])
set(gca,'Fontsize',14)
legend([h1(1) h2(1) h3(1) h4(1)], '2-point fit Diluted','2-point fit Undiluted','Linear Regression fit Diluted','Linear Regression Undiluted')
xlabel('Dilution Series Exp Growth Rate','Fontsize',14)
ylabel('Model Growth Rate','Fontsize',14)
