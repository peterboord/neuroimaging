function [fixedPeakTimes,fixedPeakLocs]=fixPhysioTimes(peakTimes,normPhysio,peakLocs,sampleTimes)

dbstop if error
if nargin>1
    peakTimes=peakTimes(:);
    normPhysio=normPhysio(:);
    peakLocs=peakLocs(:);
    sampleTimes=sampleTimes(:);
end

physioFs=(peakLocs(end)-peakLocs(1))/(peakTimes(end)-peakTimes(1));
disp(['physioFs is ',num2str(physioFs),' Hz']);
curMaxDiff=Inf;
maxDiff=0;
fixedPeakTimes=peakTimes;
% figure
noError=true;
while maxDiff<curMaxDiff && noError
%     subplot(3,1,1)
%     plot(peakTimes(1:end-1),diff(peakTimes));
%     subplot(3,1,2);
%     fixedPeakLocs=round(fixedPeakTimes*physioFs+1);
%     plot(sampleTimes,normPhysio,'b',fixedPeakTimes,normPhysio(fixedPeakLocs),'r*',peakTimes,normPhysio(peakLocs),'g*');
%     subplot(3,1,3)
%     plot(fixedPeakTimes(1:end-1),diff(fixedPeakTimes));

    diffAll=diff(fixedPeakTimes);
    medianDiffAll=median(diffAll);
    disp(['median difference is ',num2str(medianDiffAll),' sec, i.e. ',num2str(60./medianDiffAll),' bpm']);
    [maxDiff,outlierIdx]=max(abs(diffAll-medianDiffAll));
    disp(['maxDiff is ',num2str(maxDiff),' sec']);
    if maxDiff<curMaxDiff
        curMaxDiff=maxDiff;
    end
    if outlierIdx>=numel(diffAll)
        noError=false;
    else
        if diffAll(outlierIdx+1)<0.5*medianDiffAll
            % false peak - delete peak farthest from mid-point
            midPoint=mean([fixedPeakTimes(outlierIdx-1),fixedPeakTimes(outlierIdx+2)]);
            [~,peakToDeleteIdx]=max(abs(midPoint-fixedPeakTimes(outlierIdx:outlierIdx+1)));
            fixedPeakTimes(outlierIdx+peakToDeleteIdx-1)=[];
%         elseif diffAll(outlierIdx)>1.5*medianDiffAll
%             % missing peak(s) - add evenly spaced peak(s)
%             nrPeaksMissing=round((fixedPeakTimes(outlierIdx+1)-fixedPeakTimes(outlierIdx))/medianDiffAll)-1;
%             missingIsi=(fixedPeakTimes(outlierIdx+1)-fixedPeakTimes(outlierIdx))/(nrPeaksMissing+1);
%             fixedPeakTimes=[fixedPeakTimes(1:outlierIdx);...
%                 fixedPeakTimes(outlierIdx)+missingIsi*(1:nrPeaksMissing)';...
%                 fixedPeakTimes(outlierIdx+1:end)];
        end
    end
end
fixedPeakLocs=round(fixedPeakTimes*physioFs+1);
end