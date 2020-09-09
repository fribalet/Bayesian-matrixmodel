%Confidence Intervals on Slope for growth rate estimation:

%the formula for calculating the CI is b1 +/- t(0.025)SE(b1)

%the standard error of the slope is:
% [m n]=size(Y);
% SS=sum(sum((Y-repmat((dilution(1:2)*b3(2) + b3(1)),3,1)).^2)); %squared residuals from fitted line
% SSstd=sqrt(SS/(m*n -2));
% xbar=(m*dilution(1)+m*dilution(2))/(n*m);
% xSS=sum(sum((repmat([dilution(1) dilution(2)],3,1)-xbar).^2));
% SE=SSstd./sqrt(xSS);
% 
% upperCI=b1(3) + tinv(0.975,4)*SE;
% lowerCI=b1(3) - tinv(0.975,4)*SE;


%the formula for calculating the CI for the intercept b0 is:
%b0 +/- t(0.025, n-2)s(b0), where s(b0) = sqrt(MSE*(1/n + xbar2/(sum(xi-xbar)2)

%for 2pt method:
clear m n SSE MSE xbar xSS sb0 

DS_T2=DS11_T48;
DS_T1=DS11_T24;
dilution=dilutions;
% Y=log((DS_T2)./(DS_T1));
%%
%do the resgression:
Y=log((DS_T2(:,1:2))./(DS_T1(:,1:2)));
%exclude an outlier point for DS10 day1:
% Y(3,1)=NaN;
x=repmat(dilution(1:2),3,1);
x=x(:);
X=[ones(length(x),1) x]; %add coefficent to do regression

[b2pt,~,~,~,stats] = regress(Y(:),X)



[m n]=size(Y);
SSE=sum(sum((Y-repmat((dilution(1:2)*b2pt(2) + b2pt(1)),3,1)).^2));
MSE=SSE/(m*n-2);
xbar=(m*dilution(1)+m*dilution(2))/(n*m);
 xSS=sum(sum((repmat([dilution(1) dilution(2)],3,1)-xbar).^2));

sb0=sqrt(MSE*((1/(m*n)) + xbar^2./xSS));
CI = tinv(0.975,(m*n)-2)*sb0

%% for whole linear regression:
%for 2pt method:
clear m n SSE MSE xbar xSS sb0 b2pt

DS_T1=DS11_T0;
DS_T2=DS11_T24;
dilution=dilutions;
%straight linear fit:
Y=log((DS_T2)./(DS_T1));
% exclude an outlier point for DS10 day1:
% Y(3,1)=NaN;
%%
x=repmat(dilution,3,1);
x=x(:);
X=[ones(length(x),1) x]; %add coefficent to do regression

[blin,~,~,~,stats] = regress(Y(:),X);

[m n]=size(Y);
SSE=nansum(nansum((Y-repmat((dilution*blin(2) + blin(1)),3,1)).^2));
MSE=SSE/(m*n-2);
xbar=sum((m*dilution))/(n*m);
 xSS=sum(sum((repmat(dilution,3,1)-xbar).^2));

sb0=sqrt(MSE*((1/(m*n)) + xbar^2./xSS));
CI = tinv(0.975,(m*n)-2)*sb0



