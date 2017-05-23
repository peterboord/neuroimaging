function plotHemoResponse(lowHi)

switch lowHi
    case 'lo'
load('/NAS_II/Home/pboord/PIC/Udall/reactionTime/hemodynamicResponse/hrLo.mat');

for roiIdx=1:60,figure('WindowStyle','Docked'),plot([mean(hrfHi(:,roiIdx,1:25),3),mean(hrfHi(:,roiIdx,26:46),3)]);end

    case 'hi'
        
    otherwise
        error('unknown case');
end