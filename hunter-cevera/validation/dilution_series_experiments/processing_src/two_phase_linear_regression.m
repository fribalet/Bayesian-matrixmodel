%2 PHASE LINEAR REGRESSION: FITTING AND TESTING

% Here we fit the model:

%Y= Bo + B1*X + B2*(X-Xc)*I(x) +e, where Xc is the change point at which the
%lines switch. I(x) is an indicator function that is 0 if X < Xc and 1 if X
%> Xc. We fit this as a 2 regressor model, with X as one regessor and
%[(X-Xc)*I(x) as the other. This will give a residual sum of squares for
%aspecified Xc. To find this Xc, we loop over all possible values and
%identify the Xc that gives the lowest RSS. Here Bo is our desired estimate
%of the intercept, and for the dilution series, the estimate of the true
%growth rate.
DS_T2=DS9_T24;
DS_T1=DS9_T0;

%  dilution=[0.2 0.4 0.6 1];
 dilution=[0.1 0.3 0.5 0.7 1];
Y=log((DS_T2)./(DS_T1)); Y=Y(:);
% Y(3,1)=NaN; %only for day 10

x=repmat(dilution,3,1);
x=x(:);

%for DS 10, day 1:
% ii=find(~isnan(Y));
% Y=Y(ii); x=x(ii);
%%
j=0; clear rec
for Xc=0.58:0.001:0.6
    
    j=j+1;
    
    indc=x;
    indc(x>=Xc)=1; indc(x<Xc)=0;
    designX=[ones(length(x),1) x (x-Xc).*indc];
    [b,bint,r] = regress(Y,designX);
    rss=sum(r.^2);
    
    rec(j,:)=[Xc b(1) b(2) b(3) rss];
end

[mm ii]=min(rec(:,end))
rss=rec(ii,end);
%% and now just the linear model:

[blin, bint, r]=regress(Y, [ones(length(x),1) x]);
rss_lin=sum(r.^2);

%% HYPOTHESIS TESTING OF ONE PHASE VS TWO PHASE:
%The statistic for testing the one-phase model against the two-phase model is:
%F = ((RSSo - RSS1)/3) / (RSS1/(n-4))
%where RSSo is the residual sum of squares for the one-phase model, 
%RSS1 is the residual sum of squares for the two-phase model, and n is the total number of observations 
%(remember:  you're using all the data and not averaging replicates).  
%Under the null hypothesis, F has an approximate F distribution with 3 degrees of freedom in the numerator and n-4.
%So you would reject the one-phase model in favor of the two-phase model 
%at significance level alpha if the value of F is greater than the upper 0.05 quantile of this F distribution.  

F_stat=((rss_lin-rss)/3)/(rss/(length(Y)-4));
%value of F distribution with 
f_crit=finv(0.95,3,length(x)-4);

p_value=fcdf(F_stat,3,length(x)-4)
%% some plots to see how the fit looks:
  clf
%
Xc=rec(ii,1); b=rec(ii,2:4);
ii=find(dilution < Xc);
jj=find(dilution >= Xc);

plot(dilution,log((DS_T2)./(DS_T1)),'k*')
hold on
plot([dilution(ii) Xc],b(1)+b(2)*[dilution(ii) Xc],'c--')
plot([Xc dilution(jj)],b(1)+b(2)*[Xc dilution(jj)]+b(3)*([Xc dilution(jj)]-Xc),'g--')


%% with intercept:
clf
ii=find(dilution < Xc);
jj=find(dilution >= Xc);

plot(dilution,log((DS_T2)./(DS_T1)),'k*')
hold on
plot([0 dilution(ii) Xc],b(1)+b(2)*[0 dilution(ii) Xc],'c--')
plot([Xc dilution(jj)],b(1)+b(2)*[Xc dilution(jj)]+b(3)*([Xc dilution(jj)]-Xc),'g--')



%% and now a search for the CIs for the intercept

%do a search over different b0 parameters, find logL that corresponds to
%desired X2 significance level:


%LR=2*(logL0 - logL1)=n*(log(RSS_n0)-log(RSS_n1)) ~ X2 distribution with 1
%degree of freedom?
%%
% RSS_n1=((Y-designX*b)'*(Y-designX*b))./length(Y);
% logL1=-0.5*length(Y)*log(RSS_n1);
% LR=2*(logL0-logL1);
[b_hat, bint, r]=regress(Y, [ones(length(x),1) x]);
RSS_n0=sum(r.^2)./length(Y)

% RSS_n0=((Y-designX*b_hat)'*(Y-designX*b_hat))./length(Y);
logL0=-0.5*length(Y)*log(RSS_n0)
totest=logL0-chi2inv(0.95,1)


%% to find CIs using LR for single regression:
designX=x;

rec=[];
for j=b_hat(1):0.001:1.5*b_hat(1)
    b=j;
    b_star=inv(designX'*designX)*designX'*(Y-b); %optimized slope given fixed intercept
    RSS_n1=(((Y-b)-b_star.*designX)'*((Y-b)-b_star.*designX))./length(Y);
    logL1=-0.5*length(Y)*log(RSS_n1);
    rec=[rec; logL1];
    
end
   
%%
tt=b_hat(1):0.001:1.5*b_hat(1);
ii=find(rec <= totest);

CI=tt(ii(1));

%%
% for one _phase linear regression:
%X2 titles: DS exp#, day in exp, intercept, slope1, slope2, change point, RSS, CI one-sided; 
X2_CI_record=[9 1 b_hat(1) b_hat(2) NaN NaN length(Y)*RSS_n0 CI-b_hat(1); X2_CI_record]

%% for the two-phase linear regression:

Xc=ds_twophase_rec(1,3);
indc=x;
indc(x>=Xc)=1; indc(x<Xc)=0;
designX=[ones(length(x),1) x (x-Xc).*indc];
[b_hat,bint,r] = regress(Y,designX);
RSS_n0=sum(r.^2)./length(Y)

% RSS_n0=((Y-designX*b_hat)'*(Y-designX*b_hat))./length(Y);
logL0=-0.5*length(Y)*log(RSS_n0)
totest=logL0-chi2inv(0.95,1)

%%
designX=[x (x-Xc).*indc];

rec=[];
for j=b_hat(1):0.001:1.5*b_hat(1)
    b=j;
    b_star=inv(designX'*designX)*designX'*(Y-b); %optimized slope given fixed intercept
    RSS_n1=(((Y-b)-designX*b_star)'*((Y-b)-designX*b_star))./length(Y);
    logL1=-0.5*length(Y)*log(RSS_n1);
    rec=[rec; logL1];
    
end

%%
tt=b_hat(1):0.001:1.5*b_hat(1);
ii=find(rec <= totest);

CI=tt(ii(1));
%%
X2_CI_record=[6 2 b_hat(1) b_hat(2) b_hat(3) Xc length(Y)*RSS_n0 CI-b_hat(1); X2_CI_record]