function y=smooth2(a,x);
% SMOOTH y=smooth(a,x)
% Smoothing of data in columns a over every x data points
asize=size(a,1);
for rep=1:asize
		if rep<x/2+1	
			y(rep,:)=nanmean(a(1:x,:));
		elseif rep>asize-x/2-1
			y(rep,:)=nanmean(a(asize-x:asize,:));
		else
			y(rep,:)=nanmean(a(rep-x/2:rep+x/2,:));
		end
end
