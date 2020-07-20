function [prob]=loglike_DMN_14params_phours_plateau(Einterp,N_dist,theta,volbins,hr1,hr2,dt,obstimes,save_flag,savepath)

%Calculates negative log likelihood of a set of parameters given a day of
%observations (counts of cells in each size class for each hour). Model version is a two subpopulation model structure
% as described in Hunter-Cevera et al. 2014. Subpopulations are
% distinguished based on starting mean volume. The subpopulation with the
% smaller volume is referred to subpopn 1

%The inputs are as follows:

%Einterp - interpolated light data for every 10 min of the day (W/m2)
%N_dist - number of counts of cells in each size class as specified by volbins
%volbins - cell size classes (micrometers cubed)
%hr1 and hr2 refer to the starting and ending hour of the portion of day
%you want to fit. In the paper, we've used hr1=7 hours after dawn and
%hr2=25 (run till end of day).
%theta - set of parameters, described below:

gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1
b1=theta(2);  %shape parameter division function, subpopn 1
E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1
dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1
gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2
b2=theta(6); %shape parameter division function, subpopn 2
E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2
dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2
f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1
m1=theta(10); %mean volume for starting cell size distribution, subpopn 1
m2=theta(11); %mean volume for starting cell size distribution, subpopn 2
sigma1=theta(12); %variance parameter for starting cell size distributions for popn 1
sigma2=theta(13); %variance parameter for starting cell size distributions for popn 2
s=theta(14); %overdispersion parameter for the Dirichlet-multinomial distribution


gmax1=theta(1);
b1=theta(2);
E_star1=theta(3);
dmax1=theta(4);
gmax2=theta(5);
b2=theta(6);
E_star2=theta(7);
dmax2=theta(8);
f=theta(9); %fraction of starting distribution
m1=theta(10);
m2=theta(11);
sigma1=theta(12);
sigma2=theta(13);
s=100*theta(14); %scaled to help solver shrink distance between small numbers and very large numbers encountered in theta vector

q=hr2-hr1;
m=length(volbins);
ts=hr1-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create all B matrices for the hours:

B_day1=zeros(m,m,q);
B_day2=zeros(m,m,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1,dt,ts);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2,dt,ts);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y1=normpdf(1:57,m1,sigma1);
% y2=normpdf(1:57,m2,sigma2);
y1=normpdf(1:m,m1,sigma1);
y2=normpdf(1:m,m2,sigma2);

startindex = find(hr1 == obstimes);
if length(startindex) == 0
    error('Invalid value of hr1: no observation with this time')
end

Nt1=f*sum(N_dist(:,startindex))*(y1./sum(y1));
Nt2=(1-f)*sum(N_dist(:,startindex))*(y2./sum(y2));

Nt1=Nt1';
Nt2=Nt2';

simdist(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour has no zeros!

for t=1:q

    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           %project forward with the numbers
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);

    simdist(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood

    if any(isnan(simdist(:,t+1))) %just in case
        disp(['DMN 14param...simdist has a nan? theta:' num2str(theta)])
%         keyboard
    end

end

if save_flag
    file = strcat(savepath, 'fitted2.csv');
    writematrix(simdist, file);
end
%% Now calculate the log likelihood using the Dirichlet Multinomial distribution: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specifiy the expected distribution:
alpha=s*simdist;
TotN=sum(N_dist(:,startindex:end));

fitlength = length(obstimes)-startindex+1;
logL=zeros(fitlength,1);
for t=1:fitlength
    C = gammaln(s) - gammaln(TotN(:,t)+s); %constant out in front
    logL(t)=C+sum(gammaln(N_dist(:,t+startindex-1)+alpha(:,t)) - gammaln(alpha(:,t)));
end

% if any(isnan(logL))
%     keyboard
% end

prob=-sum(logL); %negative for the minimiation routine fmincon
