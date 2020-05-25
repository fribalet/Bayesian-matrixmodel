function [plot_handles]=mvco_plot(x,y,marker,color,ln_width,current_axis,gapsize)
%[plot_handles]=mvco_plot(x,y,marker,color,ln_width,current_axis,gapsize)
%plot script that breaks up the data into chunks separted by a specified
%distance (gapsize), such that if want to plot data as lines, the data gaps
%actaully appear as gaps:

temp=abs(diff(x)); %this is usually time
jj=find(temp > gapsize); %find gaps that are larger than 1 day
%setup the indexing for thedataslices:
ii=nan(2*length(jj),1);
ii(1:2:end-1)=jj;
ii(2:2:end)=jj+1;
ii=[1; ii; length(y)];

plot_handles=[];

for k=1:2:(length(ii)-1)
    yslice=y(ii(k):ii(k+1));
    xslice=x(ii(k):ii(k+1));
    p=plot(current_axis,xslice,yslice,marker,'color',color,'linewidth',ln_width);
    plot_handles=[plot_handles;p];
end
