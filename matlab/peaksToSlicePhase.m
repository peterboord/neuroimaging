function slicePhases=peaksToSlicePhase(peakTimes,tr,nslice,nvol)

slicePhases=zeros(nslice,nvol);
% augment with pre- and post-beat estimate
nrToAdd=ceil(peakTimes(1)/(peakTimes(2)-peakTimes(1)));
for addNr=1:nrToAdd
    peakTimes=[2*peakTimes(1)-peakTimes(2);peakTimes]; %#ok<AGROW>
end
nrToAdd=ceil((nvol*tr-peakTimes(end))/(peakTimes(end)-peakTimes(end-1)));
for addNr=1:nrToAdd
    peakTimes=[peakTimes;2*peakTimes(end)-peakTimes(end-1)]; %#ok<AGROW>
end
for volNr=1:nvol
    for sliceNr=1:nslice
        t=((volNr-1)*nslice+sliceNr-1)*tr/nslice;
        peakTimeBefore=peakTimes(find(peakTimes<=t,1,'last'));
        peakTimeAfter=peakTimes(find(peakTimes>t,1,'first'));
        slicePhases(sliceNr,volNr)=2*pi*(t-peakTimeBefore)/(peakTimeAfter-peakTimeBefore);
    end
end