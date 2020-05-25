%% create a structure movied that for each day that shows:
%observed data,
% best model fit and parameter values
%along with likelihood values of all modelruns:

%make a matlab movie structure or avi movie of the model fits:

% clear all
% close all



% for year2do=2008:2009;


% %avi file:
% eval(['writerobj1 = VideoWriter(''modelfits_' num2str(year2do) '.avi'');'])
% open(writerobj1);
%
%     switch year2do
%         case 2003
%             filelabel='May';
%         case 2004
%             filelabel='Apr';
%         case 2005
%             filelabel='Apr';
%         case 2006
%             filelabel='May';
%         case 2007
%             filelabel='Mar';
%         otherwise
%             filelabel='Jan';
%     end
%
disp(['Making modelfit movie for: ' num2str(year2do)])

filelist = dir([setupdays_path 'day*data.mat']);

load([modelres_path 'mvco_14par_dmn_' num2str(year2do) '.mat']);
eval('MR=modelresults;')
eval('allMR=allmodelruns;')

%for each day:
clf
set(gcf,'Position',[119         257        1705         673])
hr1=7; hr2=25; xsp=0.035; ysp=0.1;
%%
for j=1:length(filelist)
    
    filename=filelist(j).name;
    day=str2num(filename(4:9));
    eval(['load ' setupdays_path filename])
    
    jj=find(MR(:,1)==day);
    %
    if ~isempty(jj)
        
        theta=MR(jj,2:15);
        
        %Fix and Interpolate Light Data:
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(find(isnan(Einterp))) = 0;
        
        if size(N_dist,2) < 25
            m=size(N_dist,2);
            N_dist=[nan(57,25-m) N_dist];
            Vhists=[nan(57,25-m) Vhists];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %project forward plus growth and division functions:
        
        [dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample_plt(Einterp,N_dist,theta,volbins,hr1,hr2);
        dirsample=dirsample./repmat(sum(dirsample,1),57,1);
        dirsample=[nan(57,6) dirsample];
        simdist=[nan(57,6) simdist];
        
        del1=(theta(4).*volbins.^theta(2))./(1+(volbins.^theta(2)));
        y1=theta(1)*ones(size(Einterp));
        ind=find(Einterp < theta(3));
        y1(ind)=(theta(1)/theta(3)) * Einterp(ind);
        
        del2=(theta(8).*volbins.^theta(6))./(1+(volbins.^theta(6)));
        y2=theta(5)*ones(size(Einterp));
        ind=find(Einterp < theta(7));
        y2(ind)=(theta(5)/theta(7)) * Einterp(ind);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %for each day in the year2do, see how model distributions compare
        %with observations:
        
        % avi file with no output ot screen:
        %             fig2=figure(2);
        %             clf(fig2)
        %             set(gcf,'Position',[ 93 274  1649 655],'visible','off')
        subplot(2,7,7,'replace')
        subplot_tight(2,7,7,[ysp xsp-0.01])
        title('Parameter Values','fontsize',14)
        text(0.0,0.6,{['gmax 1: ' num2str(theta(1))];...
            ['b 1: ' num2str(theta(2))];
            ['E* 1: ' num2str(theta(3))];
            ['dmax 1: ' num2str(theta(4))];
            ['vol 1: ' num2str(theta(10))];
            ['sigma1: ' num2str(theta(12))];
            ['proportion: ' num2str(theta(9))];
            ['mu2: ' num2str(MR(jj,18))];
            ['ending p: ' num2str(MR(jj,20))]})
        
        
        text(0.6,0.6,{['gmax 2: ' num2str(theta(5))];...
            ['b 2: ' num2str(theta(6))];
            ['E* 2: ' num2str(theta(7))];
            ['dmax 2: ' num2str(theta(8))];
            ['vol 2: ' num2str(theta(11))];
            ['sigma2: ' num2str(theta(13))];
            ['proportion: ' num2str(1-theta(9))];
            ['mu2: ' num2str(MR(jj,19))];
            ['ending p: ' num2str(MR(jj,21))]})
        
        text(0.0, 0.1,{['s: ' num2str(100*theta(14))];
            ['mu: ' num2str(MR(jj,17))]});
        set(gca,'visible','off')
        
        subplot(2,7,1,'replace')
        subplot_tight(2,7,1,[ysp xsp])
        plot(Edata(:,1),Edata(:,2),'.-'), hold on
        plot(time,Einterp,'k.')
        xlabel('Time')
        ylabel('Edata')
        
        subplot(2,7,2,'replace')
        subplot_tight(2,7,2,[ysp xsp])
        h0=imagesc(1:25,1:57,Vhists); set(gca,'Ydir','normal')
        set(h0,'AlphaData',~isnan(Vhists));
        xlabel('Time')
        ylabel('Cell size class')
        title('Observed Data')
        title([datestr(day) ' : ' num2str(day)])
        
        subplot(2,7,3,'replace')
        subplot_tight(2,7,3,[ysp xsp])
        h=imagesc(1:25,1:57,simdist);
        set(h,'AlphaData',~isnan(simdist)), set(gca,'Ydir','normal')
        xlabel('Time')
        ylabel('Cell size class')
        title('MLE Model Fit - no Dir. sample')
        
        subplot(2,7,4,'replace')
        subplot_tight(2,7,4,[ysp xsp])
        h1=imagesc(1:25,1:57,dirsample);
        set(h1,'AlphaData',~isnan(dirsample)), set(gca,'Ydir','normal')
        xlabel('Time')
        ylabel('Cell size class')
        title('Sample from Dirichlet')
        
        subplot(2,7,5,'replace'), hold on
        subplot_tight(2,7,5,[ysp xsp]), hold on
        [~,iy]=sort(y1);
        plot(Einterp(iy),y1(iy),'.-','color',[0 0.5 1],'markersize',8);
        [~,iy]=sort(y2);  set(gca,'box','on')
        plot(Einterp(iy),y2(iy),'.-','color',[0 0 0.8],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Fraction growing cells')
        xlabel('Edata')
        
        subplot(2,7,6,'replace')
        subplot_tight(2,7,6,[ysp xsp]), hold on
        temp=allMR{jj};
        [~, ib]= sort(temp(:,15));
        plot(1:size(temp,1),temp(ib,16),'-','color',[0.6 0.6 0.6]), hold on
        plot(1:size(temp,1),temp(ib,16),'k.')
        xlim([0 size(temp,1)])
        YL=get(gca,'ylim');  ylim([0 YL(2)+0.1]);  set(gca,'box','on')
        line([5 5],[0 YL(2)+0.1],'color',[0.6 0.6 0.6])
        ylabel('Division rate')
        xlabel('# model runs')
        
        subplot(2,7,8,'replace')
        subplot_tight(2,7,8,[ysp xsp])
        h1=plot(1:57,Vhists(:,hr1:hr1+4),'color',[0.5 0.5 0.5]); hold on
        h2=plot(1:57,simdist(:,1:5),'color',[0 0 0]);
        h3=plot(1:57,Vt1(:,1:5),'color',[0 0.5 1]);
        h4=plot(1:57,Vt2(:,1:5),'color',[0 0 0.8]);
        plot(1:57,Vt1(:,1),'color',[1 0.5 0]);
        plot(1:57,Vt2(:,1),'color',[1 0 0]);
        title('Hours 7-11')
        legend([h1(1); h2(1); h3(1); h4(1)],'Obs','Model','popn1','popn2','location','NorthWest')
        xlabel('Size class')
        ylabel('Proportion')
        
        subplot(2,7,9,'replace')
        subplot_tight(2,7,9,[ysp xsp])
        plot(1:57,Vhists(:,hr1+5:hr1+9),'color',[0.5 0.5 0.5]), hold on
        plot(1:57,simdist(:,6:10),'color',[0 0 0])
        plot(1:57,Vt1(:,6:10),'color',[0 0.5 1])
        plot(1:57,Vt2(:,6:10),'color',[0 0 0.8])
        title('Hours 12-16')
        xlabel('Size class')
        ylabel('Proportion')
        
        subplot(2,7,10,'replace')
        subplot_tight(2,7,10,[ysp xsp])
        plot(1:57,Vhists(:,hr1+10:hr1+14),'color',[0.5 0.5 0.5]), hold on
        plot(1:57,simdist(:,11:15),'color',[0 0 0])
        plot(1:57,Vt1(:,11:15),'color',[0 0.5 1])
        plot(1:57,Vt2(:,11:15),'color',[0 0 0.8])
        title('Hours 17-21')
        xlabel('Size class')
        ylabel('Proportion')
        
        subplot(2,7,11,'replace')
        subplot_tight(2,7,11,[ysp xsp])
        plot(1:57,Vhists(:,hr1+15:hr1+18),'color',[0.5 0.5 0.5]), hold on
        plot(1:57,simdist(:,16:19),'color',[0 0 0])
        plot(1:57,Vt1(:,16:19),'color',[0 0.5 1])
        plot(1:57,Vt2(:,16:19),'color',[0 0 0.8])
        title('Hours 22-25')
        xlabel('Size class')
        ylabel('Proportion')
        
        subplot(2,7,12,'replace'), hold on
        subplot_tight(2,7,12,[ysp xsp]), hold on
        plot(volbins,del1,'.-','color',[0 0.5 1],'markersize',8);
        plot(volbins,del2,'.-','color',[0 0 0.8],'markersize',8);
        ylabel('Fraction dividing cells')
        xlabel('Size bins')
        
                
        subplot(2,7,14,'replace')
        subplot_tight(2,7,14,[ysp xsp]), hold on
        plot(MR(:,1),MR(:,17),'.')
        plot(MR(jj,1),MR(jj,17),'rp')
        datetick('x','mm/dd')
        ylabel('Division rate')
        set(gca,'box','on')
        
        subplot(2,7,13,'replace')
        subplot_tight(2,7,13,[ysp xsp]), hold on
        temp=allMR{jj};
        [~, ib]= sort(temp(:,15));
        plot(1:size(temp,1),temp(ib,15),'-','color',[0.6 0.6 0.6]), hold on
        plot(1:size(temp,1),temp(ib,15),'.')
        line([5 5],ylim,'color',[0.6 0.6 0.6])
        ylabel('-log L')
        xlabel('# model runs')
        set(gca,'box','on')

        
       %% 
        F1(j)=getframe(gcf); %for visible output
        %F1(w) = im2frame(zbuffer_cdata(fig1)); %if don't want figure to appear on screen
        
        %             %for avi file:
        %             F2 = im2frame(zbuffer_cdata(fig1));
        %             writeVideo(writerobj1,F2);
    else
        clf, text(0.5,0.5,['didn''t make day: ' datestr(day) ' : ' num2str(day) '?'],'fontsize',20)
        set(gca,'visible','off')
        F1(j)=getframe(gcf); %for visible output
    end
end

eval(['modelfits' num2str(year2do) '=F1;'])
eval(['save ' modelres_path 'modelfits' num2str(year2do) ' modelfits' num2str(year2do) '' ])
%
disp(['Finished ' num2str(year2do) '!'])
%   %avi file:
%   close(writerobj1);

%    clearvars -except hr1 hr2 year2do
% end
