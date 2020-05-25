function [dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample_plt(Einterp,N_dist,theta,volbins,hr1,hr2)

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
s=100*theta(14);

q=hr2-hr1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create all B matrices for the hours:

B_day1=zeros(57,57,q);
B_day2=zeros(57,57,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1=normpdf(1:57,m1,sigma1);
y2=normpdf(1:57,m2,sigma2);

Nt1=f*sum(N_dist(:,hr1))*(y1./sum(y1));
Nt2=(1-f)*sum(N_dist(:,hr1))*(y2./sum(y2));

Nt1=Nt1';
Nt2=Nt2';
Vt1=Nt1./sum(Nt1+Nt2);
Vt2=Nt2./sum(Nt1+Nt2);
simdist(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2);  %Only if starting hour has no zeros!

for t=1:q

    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           %project forward with the numbers
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    Vt1(:,t+1)=Nt1(:,t+1)./sum(Nt1(:,t+1)+Nt2(:,t+1));
    Vt2(:,t+1)=Nt2(:,t+1)./sum(Nt1(:,t+1)+Nt2(:,t+1));
    simdist(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood

    if any(isnan(simdist(:,t+1))) %just in case
        disp(['DMN 13param...simdist has a nan? theta:' num2str(theta)])
%         keyboard
    end

end

%simdist is then used to generate aample from the Dirichlet distribution:

dirsample=zeros(57,q);
for i=1:q+1
dirsample(:,i)=gamrnd(s*simdist(:,i),1);
dirsample(:,i)=dirsample(:,i)./sum(dirsample(:,i));
dirsample(:,i)=mnrnd(sum(N_dist(:,hr1-1+i)),dirsample(:,i));
end

end
