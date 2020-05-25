%setup days quality control on how days were constructed:
MT={'2003' 'May';
    '2004' 'Apr';
    '2005' 'Apr';
    '2006' 'May';
    '2007' 'Mar';
    '2008' 'Jan';
    '2009' 'Jan';
    '2010' 'Jan';
    '2011' 'Jan';
    '2012' 'Jan';
    '2013' 'Jan';
    '2014' 'Jan';
    '2015' 'Jan';
    '2016' 'Jan'};

monthtag=MT{year2do-2002,2};

if ~isempty(strfind(computer,'MC'))
        daypath=['/Volumes/Lab_data/MVCO/FCB/MVCO_' monthtag num2str(year2do) '/model/input_beadmean_July2016/'];
        savepath=daypath;
else
         daypath=['\\sosiknas1\lab_data\MVCO\FCB\MVCO_' monthtag num2str(year2do) '\model\input_beadmean_July2016\'];
        savepath=daypath;    
end

filelist=dir([daypath 'day*data.mat']);

%% in a simple matlab movie form (results in a structure):
figure, set(gcf,'Position',  [ 82         428        1524         521])
count=0;

for j=1:length(filelist)
    
    count=count+1;
    
    filename=filelist(j).name;
    eval(['load ' daypath filename])
    
    clf
    subplot(1,3,2,'replace')
    imagesc([1 25],[1 57],Vhists)
    set(gca,'YDir','normal');
    %pcolor(Vhists), hading flat
    caxis([0 0.10])
    colorbar
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Cell Volume','fontsize',14)
    title({['Day: ' filename(4:9) ' : ' datestr(str2num(filename(4:9)))]; 'Vhists'},'fontsize',14)
    
    subplot(133)
    imagesc([1 25],[1 57],N_dist)
    set(gca,'YDir','normal');
    %         pcolor(N_dist)
    %         shading flat
    colorbar %caxis([0 10000])
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Cell Volume','fontsize',14)
    t1=title('N_dist','fontsize',14);
    set(t1,'interpreter','none')
    
    subplot(131)
    plot(Edata(:,1), Edata(:,2),'.-','color',[0 0.3 0.8],'markersize',10,'linewidth',1.5)
    xlabel('Hours After Dawn','fontsize',14)
    ylabel('Radiation','fontsize',14)
    ylim([0 1000])
    
    F1(count) = getframe(gcf);
end

notes='to play, type implay(syn_dist20XX) in MATLAB window';
eval(['syn_dist' num2str(year2do) '=F1;'])
eval(['save ' savepath 'setup_days' num2str(year2do) '.mat syn_dist' num2str(year2do) ' notes'])

close all
clear F1
