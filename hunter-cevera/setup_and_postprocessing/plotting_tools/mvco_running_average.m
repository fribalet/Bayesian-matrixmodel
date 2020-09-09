function totalrunavg=mvco_running_average(time,data,windowsize,gapsize)
%totalrunavg=mvco_running_average(time,data,windowsize,gapsize)
%function to calculate a running average of abundance data
%note that this function is written to handle MVCO data, such that breaks
%in the dataset are dealt with and do not carry over into the averages of
%different deployments for example
%gapsize is the time interval at which to start considering different
%blocks of data as different chunks (needs to have units of days!)


if floor(windowsize/2) ~= windowsize/2 %means windowsize is odd
    windowsize=windowsize-1; %needs to be even for indexing below
end

temp=abs(diff(time)); %this is usually the cellresultsall
jj=find(temp > gapsize); %find gaps that are larger than 1 day
%setup the indexing for thedataslices:
ii=nan(2*length(jj),1);
ii(1:2:end-1)=jj;
ii(2:2:end)=jj+1;
ii=[1; ii; length(data)];
totalrunavg=nan(length(data),1);

for k=1:2:(length(ii)-1)
    dataslice=data(ii(k):ii(k+1));
    runavg=nan(length(dataslice),1);
    
    if length(dataslice) >= windowsize
        for q=1:length(dataslice)
            if q <= windowsize/2
                runavg(q)=nanmean(dataslice(1:q+windowsize/2));
            elseif q >= length(dataslice)-windowsize/2
                runavg(q)=nanmean(dataslice(q-windowsize/2:end));
            else
                runavg(q)=nanmean(dataslice(q-windowsize/2:q+windowsize/2));
            end
        end
    elseif length(dataslice) <= windowsize %adaptively reduce window size for this section, based on the size of the section:
        %keyboard %make sure that this new window size is okay!
        for q=1:length(dataslice)
            
            if floor(length(dataslice)/2) ~= length(dataslice)/2 %means dataslice is odd length
                windowsizeb=(length(dataslice)+1)/2; %needs to be even for indexing below
            else
                windowsizeb=length(dataslice)/2;
            end
            disp(['Dataslice length: ' num2str(length(dataslice)) ' Adaptive window: ' num2str(windowsizeb)])
            if q <= windowsizeb/2
                runavg(q)=nanmean(dataslice(1:q+windowsizeb/2));
            elseif q >= length(dataslice)-windowsizeb/2
                runavg(q)=nanmean(dataslice(q-windowsizeb/2:end));
            else
                runavg(q)=nanmean(dataslice(q-windowsizeb/2:q+windowsizeb/2));
            end
        end
    end
    
    totalrunavg(ii(k):ii(k+1))=runavg;
end