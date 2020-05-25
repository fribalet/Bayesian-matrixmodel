%exploration of grazing functional responses and dilutions:

%the grazing rate is proportional to the clearnace rate (vol filtered  per
%grazer per time) multiplied by the grazer concentration
%The clearance rate is equal to the ingestion rate( phyto consumer per
%grazer per time) divided by the phyto concentration ( C=I/P )

%Let I =f(P) (some functional response to P)
%then, r = k - g = k - CZ = k - Z(I/P) = k - Z/P*f(P)
%for a dilution,
%this becomes: r = k - XZ1 / XP1 *F(XP1) where X is the dilution factor
% and reduces to r = k - Z1/P1*f(XP1)

% Type I  -linear then abrupt saturation
phyto_t0=[100 300 500 700 1000];
grazers_t0=[10 30 50 70 100];

igest_curve=[5 5 5 5 5]; % this is a flat/saturated grazing curve
igest_curve=[1 2 3 5 5];
igest_curve=(3*phyto_t0)./(300 + phyto_t0)

mu=0.7*ones(1,5);
L=-(grazers_t0(end)./phyto_t0(end)).*igest_curve

%% Normalized functional response as calculated in Gallegos 1989:

%normalized funcitonal response (at dilution x) = (rx - k)/-g

figure
hold on

for q=3:length(ds_exp_mu)
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
        k=ds_exp_mu(q,4);
    else
        grazing_rate=ds_exp_mu(q,5)-log(nanmean((DS_T2(:,end)))./nanmean((DS_T1(:,end))));
        k=ds_exp_mu(q,5);
    end
    r=log(nanmean((DS_T2))./nanmean((DS_T1)));
    nfr=(r-k)/(-grazing_rate);
%     plot(dilutions,nfr,'s')
    plot(nanmean(DS_T1),nfr,'o')

%     record=[record; ds k grazing_rate log(nanmean((DS_T2))./nanmean((DS_T1)))];
end
 
%% for the picoeuks:
%normalized funcitonal response (at dilution x) = (rx - k)/-g

figure
% hold on

for q=3:length(ds_exp_mu_picoeuks)
    ds=ds_exp_mu_picoeuks(q,1);
    if ds_exp_mu_picoeuks(q,2)==1
        eval(['load DS' num2str(ds) '_rawdata'])
        eval(['DS_T1= DS' num2str(ds) '_picoeuks_T0;'])
        eval(['DS_T2= DS' num2str(ds) '_picoeuks_T24;'])
    else 
        eval(['DS_T1= DS' num2str(ds) '_picoeuks_T24;'])
        eval(['DS_T2= DS' num2str(ds) '_picoeuks_T48;'])
    end
    if ~isnan(ds_exp_mu_picoeuks(q,3))
        grazing_rate=ds_exp_mu_picoeuks(q,3)-log(nanmean((DS_T2(:,end)))./nanmean((DS_T1(:,end))));
        k=ds_exp_mu_picoeuks(q,3);
    else
        grazing_rate=ds_exp_mu_picoeuks(q,4)-log(nanmean((DS_T2(:,end)))./nanmean((DS_T1(:,end))));
        k=ds_exp_mu_picoeuks(q,4);
    end
    r=log(nanmean((DS_T2))./nanmean((DS_T1)));
    nfr=(r-k)/(-grazing_rate);
    plot(dilutions,nfr,'s')
    pause
%     plot(nanmean(DS_T1),nfr,'o')

%     record=[record; ds k grazing_rate log(nanmean((DS_T2))./nanmean((DS_T1)))];
 end