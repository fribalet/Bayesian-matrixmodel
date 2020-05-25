%main script to process field data with 14-parameter model assuming a dirichlet multinomial distribution
%for likelihood calculation:

%script runs batches of random startpoints with multistart as
%this seems to be faster than running fmincon one solver run at time to

% clear all
% close all

tic; % starts timer

% select dataset to use
% data = 'simulated';
data = 'zinser';

%some front matter:
% time window for the model
if strcmp(data, 'zinser')
    hr1 = 1; hr2 = 49; % for Pro there is no restriction on division
    dt = 1/3;  % value we use for the Zinser data and the Bayesian model
    obstimes = 1:2:hr2; % specify the time associated with each observation
else
    hr1 = 7; hr2 = 25; % for Syn we assume no division for 6 hours
    dt = 1/6; % typical value used in most literature
    obstimes = 1:hr2; % specify the time associated with each observation
end

save_flag = false; % whether or not to save the fitted size distribution

% parameter naming and notes
restitles={'day';'gmax1';'b1';'E*1';'dmax1';'gmax2';'b2';'E*2';'dmax2';'proportion';'m1';'m2';'sigma1';'sigma2';'s';'-logL';'mu';'mu1';'mu2';'ending proportion 1';'ending proportion 2'; 'exitflag';'number solver runs'};
notes='E* bounds are from 0 to max(Einterp)';

% optimization setup
ms=MultiStart('Display','off','TolX',1e-5,'UseParallel','always','StartPointsToRun','bounds');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 30000,'MaxFunEvals',100000);
icsTol=0.2;
tolvec=[0.01 0.01 100 0.005 0.01 0.01 100 0.005 0.01 0.5 0.5 0.5 0.5 10];

% File paths for the datasets we use
if strcmp(data, 'zinser')
    eval(['pathname=''~/git_environment/Bayesian-matrixmodel/hunter-cevera/zinser/data/'';']) % Zinser data
    eval(['savepath=''~/git_environment/Bayesian-matrixmodel/hunter-cevera/zinser/output/'';']) % Zinser model output
    filelist = dir([pathname 'Zinser_data.mat']); % Zinser data
else
    eval(['pathname=''~/git_environment/Bayesian-matrixmodel/hunter-cevera/example/data/'';']) % Simulated data
    eval(['savepath=''~/git_environment/Bayesian-matrixmodel/hunter-cevera/example/output/'';']) % Simulated model output
    filelist = dir([pathname 'simul_Syn.mat']); % simulated data
end

%save results as they are processed in a separate directory:
if ~exist(savepath,'dir'), mkdir(savepath), end %make directory if doesn't exist yet

% setup result variables:
modelresults=zeros(length(filelist),length(restitles)-1);
allmodelruns=cell(length(filelist),2);

% For 2005 - end of October onward (732607), merging problem - skipping for
% now

%%
for filenum=1:length(filelist)

    filename=filelist(filenum).name
%     day=str2num(filename(4:9));

    %     if ~isempty(find(daylist==day))
%     disp(['optimizing day: ' num2str(day) ' file#: ' num2str(filenum)])
    eval(['load ' pathname filename])

    %Fix and Interpolate Light Data:
    nnind = find(~isnan(Edata(:,2)));
    Edata=Edata(nnind,:);
    [unqE eind]=unique(Edata(:,1));
    % If PAR length matches observation length, overwrite observation
    % times, assuming we have one PAR measurement for each observation
%     if size(N_dist,2) == size(Edata,1)
%         obstimes = Edata(:,1)
%     end
    time=Edata(1,1):(dt):Edata(end,1);
    Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
    Einterp(find(isnan(Einterp))) = 0;

    lb=-[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 5 5 1 1 1e-4]; %parameter bounds
    ub=[1 15 max(Einterp) 1 1 15 max(Einterp) 1 0.5 50 50 15 15 1e4];

    a1=-1*eye(14); %set parameter bounds to be interpretted by fmincon
    a2=eye(14);
    A=zeros(28,14);
    A(1:2:27,:)=a1;
    A(2:2:28,:)=a2;

    B=zeros(27,1);
    B(1:2:27)=lb;
    B(2:2:28)=ub;

    %starting conditions:
    x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.5*rand 30*rand+20 30*rand+20 10*rand+2 10*rand+2 1e4*rand];
    %random start points:
    tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.5*rand(40,1) 30*rand(40,1)+20 30*rand(40,1)+20 10*rand(40,1)+2 10*rand(40,1)+2 1e4*rand(40,1)]);

    problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,N_dist,theta,volbins,hr1,hr2,dt,obstimes,false,savepath),'Aineq',A,'bineq',B,'options',opts);
    [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);

    %open up the soln structure:
    temp=zeros(40,17);
    start_points=zeros(40,14);
    temp=zeros(40,17);
    c=1;

    for j=1:length(soln)
        %check to see if all start points led to an individual solution or
        %not (MultiSTart will only return unique solutions)
        g=cell2mat(soln(j).X0);
        if length(g)==14 %only one start_point led to that solution
            start_points(c,:)=g;
            temp(c,1:14)=soln(j).X;
            temp(c,15)=soln(j).Fval;
            temp(c,16)=growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2,dt,obstimes);
            temp(c,17)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/14;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
            temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
            temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2,obstimes),num,1);
            temp(c:c+num-1,17)=repmat(soln(j).Exitflag,num,1);
            c=c+num;
        end
    end
    %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);

    largepopn=zeros(size(temp,1),7);
    smallpopn=zeros(size(temp,1),7);
    for h=1:size(temp,1)
        if temp(h,10) > temp(h,11)
            largepopn(h,:) = temp(h,[1:4 9 10 12]);
            smallpopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
        else %11 > 10
            largepopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
            smallpopn(h,:) = temp(h,[1:4 9 10 12]);
        end
    end

    modelfits=[smallpopn(:,1:4) largepopn(:,1:4) smallpopn(:,5) smallpopn(:,6) largepopn(:,6) smallpopn(:,7) largepopn(:,7) temp(:,14:end)];
    start_points=start_points(qq,:);
    allstarts=start_points;

    %let's now ask, in the first batch run, did the solver "converge"?
    [sortlogL ii]=sort(modelfits(:,15));
    if abs(sortlogL(min(5,size(sortlogL,1)))-sortlogL(1)) < icsTol
        flag1 = 0;
    else
        disp(num2str(sortlogL(1:min(5,size(sortlogL,1)))))
        flag1 = 1;
    end;

    partol=max(modelfits(ii(1:min(5,size(sortlogL,1))),1:14))-min(modelfits(ii(1:min(5,size(sortlogL,1))),1:14));
    if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0;
    else
        flag2 = 1;
    end

    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    k=1; %batch number
    while ((flag1 || flag2) | size(modelfits,1) <= 20) && k <= 5

        disp(['k: ' num2str(k)])
        k=k+1;
        x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.2*rand 6*rand max(Einterp)*rand 0.1*rand 0.5*rand 30*rand+20 30*rand+20 10*rand+2 10*rand+2 1e4*rand]; %random start point

        tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 0.5*rand(40,1) 30*rand(40,1)+20 30*rand(40,1)+20 10*rand(40,1)+2 10*rand(40,1)+2 1e4*rand(40,1)]);

        problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,N_dist,theta,volbins,hr1,hr2,dt,obstimes,false,savepath),'Aineq',A,'bineq',B,'options',opts);
        [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);

        %open up the soln sturcutre:
        temp=zeros(40,17);
        start_points=zeros(40,14);
        temp=zeros(40,17);
        c=1;

        for j=1:length(soln)
            %check to see if all start points led to an individual solution or
            %not (MultiSTart will only return unique solutions)
            g=cell2mat(soln(j).X0);
            if length(g)==14 %only one start_point led to that solution
                start_points(c,:)=g;
                temp(c,1:14)=soln(j).X;
                temp(c,15)=soln(j).Fval;
                temp(c,16)=growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2,dt,obstimes);
                temp(c,17)=soln(j).Exitflag;
                c=c+1;
            else
                num=length(g)/14;
                start_points(c:c+num-1,:)=squeeze(reshape(g',1,14,num))';
                temp(c:c+num-1,1:14)=repmat(soln(j).X,num,1);
                temp(c:c+num-1,15)=repmat(soln(j).Fval,num,1);
                temp(c:c+num-1,16)=repmat(growth_rate(Einterp,volbins,N_dist,temp(c,1:13),hr1,hr2,dt,obstimes),num,1);
                temp(c:c+num-1,17)=repmat(soln(j).Exitflag,num,1);
                c=c+num;
            end
        end
        %just in case have rows left as zeros
        qq=find(temp(:,1)~=0);
        temp=temp(qq,:);
        start_points=start_points(qq,:);

        largepopn=zeros(size(temp,1),7); %large population has the larger mean starting volume bin
        smallpopn=zeros(size(temp,1),7);
        for h=1:size(temp,1)
            if temp(h,10) > temp(h,11)
                largepopn(h,:) = temp(h,[1:4 9 10 12]);
                smallpopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
            else %11 > 10
                largepopn(h,:) = [temp(h,5:8) 1-temp(h,9) temp(h,[11 13])];
                smallpopn(h,:) = temp(h,[1:4 9 10 12]);
            end
        end

        modelfits=[modelfits; smallpopn(:,1:4) largepopn(:,1:4) smallpopn(:,5) smallpopn(:,6) largepopn(:,6) smallpopn(:,7) largepopn(:,7) temp(:,14:end)];
        allstarts=[allstarts; start_points];

        %okay, now see after this batch run, did the solver "converge"?
        [sortlogL ii]=sort(modelfits(:,15));

        if abs(sortlogL(5)-sortlogL(1)) < icsTol
            flag1 = 0;
        else
            disp(num2str(sortlogL(1:5))) %should be 5, but occassionally get less than 5 solver runs returned...
            flag1 = 1;
        end;

        partol=max(modelfits(ii(1:5),1:14))-min(modelfits(ii(1:5),1:14));
        if sum(abs(partol) < tolvec)==14 || sum((abs(partol./modelfits(ii(1),1:14)) < 0.05))==14 %either the modelfits are within an absolute tolerance or within a relative tolerance
            flag2 = 0;
        else
            flag2 = 1;
        end
        disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    end  %while loop

    [s jj]=sort(modelfits(:,15));
    xmin=modelfits(jj(1),1:14);
    fmin=modelfits(jj(1),15);
    exitflag=modelfits(jj(1),17);

    [mu mu1 mu2 p1 p2]=growth_rate(Einterp,volbins,N_dist,xmin(1:13),hr1,hr2,dt,obstimes);

    modelresults(filenum,:)=[xmin fmin mu mu1 mu2 p1 p2 exitflag length(modelfits)];
    resultstable = array2table(modelresults);
    resultstable.Properties.VariableNames(1:length(modelresults)) = restitles(2:end);
    allmodelruns{filenum,1}=modelfits;
    allmodelruns{filenum,2}=allstarts;

    eval(['save ' savepath 'mvco_14par_dmn_modelresults allmodelruns'])
    eval(['save ' savepath 'mvco_14par_dmn_modelpars resultstable'])

end

% if save_flag is true, we save the final fitted size distribution
if save_flag
    negloglike_calc(Einterp,N_dist,modelresults(1:14),volbins,hr1,hr2,dt,obstimes,save_flag,savepath);
end

toc; % end timer
