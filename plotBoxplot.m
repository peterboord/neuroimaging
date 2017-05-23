dbstop if error
% load fisheriris
% x = meas(51:150,3);
% y = meas(51:150,4);
% g = species(51:150);
% 
% aoctool(x,y,g)

% group2 = ismember(g,'virginica');
% group1 = ~group2;
% X = [group2, x.*group1, x.*group2];
% s = regstats(y,X);
% [p,F] = linhyptest(s.beta, s.covb, 0, [0 0 1 -1], s.tstat.dfe)

cd '/projects2/udall/pboord/analysis/tables'
rtNocue=1;
rtCenter=2;
rtSpatial=3;
rtCon=4;
rtIncon=5;
rtAlert=6;
rtOrient=7;
rtExec=8;
c1actCon=1;
c1actIncon=2;
c1actExec=3;
c2actCon=4;
c2actIncon=5;
c2actExec=6;
c3actCon=7;
c3actIncon=8;
c3actExec=9;
c4actCon=10;
c4actIncon=11;
c4actExec=12;
c1z1c2=13;% r-thalamus
c1z2c2=14;% l-mPFC
c2z2c1=15;% lr-mPFC
c2z2c2=16;% l-striatum
c2z2c3=17;% r-precuneus
c2z2c5=18;% r-striatum
c3z1c1=19;% l-cerebellum
c4z2c1=20;% lr-precuneus
c4z2c2=21;% lr-mPFC
c4z2c3=22;% lr-striatum

%%
% figure
% pdSubj=[1:20,22:25];
% hcSubj=26:46;
% vars=[10,11,4,5]; % RIPS con, RIPS incon, RFEF con, RFEF incon
% varNames={'RIPS','RFEF'};
% % vars=[1,2,4,5,7,8,10,11];
% % varNames={'PREC','RFEF','LIPS','RIPS'};
% nrExec=floor(numel(vars)/2);
% % Activation in executive task clusters during congruent and incongruent trials.
% actFcData=csvread('actFcData/actFcData.csv');
% actFcHdrs=textread('actFcData/actFcData_colHeaders.txt','%s');
% edata=actFcData(:,vars);
% edata={edata(pdSubj,:),edata(hcSubj,:)};
% mnEdata=reshape(cell2mat(cellfun(@mean,edata,'UniformOutput',0)),2,[]);
% sdEdata=reshape(cell2mat(cellfun(@std,edata,'UniformOutput',0)),2,[]);
% szEdata=reshape(cell2mat(cellfun(@length,edata,'UniformOutput',0)),2,[]);
% seEdata=sdEdata./reshape(repmat(sqrt(szEdata),1,numel(vars))',2,[]);
% mnEdata=([mnEdata(:,1:nrExec);mnEdata(:,nrExec+1:2*nrExec)])';
% seEdata=([seEdata(:,1:nrExec);seEdata(:,nrExec+1:2*nrExec)])';
% % barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
% barweb(mnEdata,seEdata,[],varNames,[],'Executive contrast cluster','Z statistic of contrast (con/incon)',[],[],{'PD Congruent','PD Incongruent','HC Congruent','HC Incongruent'},[],'plot');
% 
% keyboard

%%
age=load('age.txt','-ascii');
rtHdrs=textread('reactionTime_baseline_colHeaders.txt','%s');
rtData=csvread('reactionTime_baseline.csv');
actFcData=csvread('actFcData/actFcData.csv');
actFcHdrs=textread('actFcData/actFcData_colHeaders.txt','%s');
fcData=csvread('fcData/fcData.csv');
fcHdrs=textread('fcData/fcData_colHeaders.txt','%s');

%%
figure
pdSubj=[1:20,22:25];
hcSubj=26:46;
set(gcf, 'Position', get(0,'Screensize'));
%% subplot(2,2,1);
subplot(2,2,1);
plotBoxplot_subplot(1,actFcData,c4z2c2,c4actExec,pdSubj,'','PD');
%% subplot(2,2,2);
subplot(2,2,2);
plotBoxplot_subplot(2,actFcData,c2actExec,c4actExec,pdSubj,'','');
%% subplot(2,2,3);
subplot(2,2,3);
plotBoxplot_subplot(3,actFcData,c4z2c2,c4actExec,hcSubj,'RIPS-mPFC FC','HC','RIPS Activation');
%% subplot(2,2,4);
subplot(2,2,4);
plotBoxplot_subplot(4,actFcData,c2actExec,c4actExec,hcSubj,'RFEF Activation','');

keyboard

%%
figure
pdSubj=1:25;
hcSubj=26:46;
set(gcf, 'Position', get(0,'Screensize'));
%% subplot(2,4,1);
subplot(2,4,1);
plotBoxplot_subplot(1,actFcData,c4z2c2,c4actExec,pdSubj,'','PD');
%% subplot(2,4,2);
subplot(2,4,2);
plotBoxplot_subplot(2,actFcData,c1actExec,c4actExec,pdSubj,'','');
%% subplot(2,4,3);
subplot(2,4,3);
plotBoxplot_subplot(3,actFcData,c2actExec,c4actExec,pdSubj,'','');
%% subplot(2,4,4);
subplot(2,4,4);
plotBoxplot_subplot(4,actFcData,c3actExec,c4actExec,pdSubj,'','');
%% subplot(2,4,5);
subplot(2,4,5);
plotBoxplot_subplot(5,actFcData,c4z2c2,c4actExec,hcSubj,'RIPS-mPFC FC','HC','RIPS Activation');
%% subplot(2,4,6);
subplot(2,4,6);
plotBoxplot_subplot(6,actFcData,c1actExec,c4actExec,hcSubj,'PREC Activation','');
%% subplot(2,4,7);
subplot(2,4,7);
plotBoxplot_subplot(7,actFcData,c2actExec,c4actExec,hcSubj,'RFEF Activation','');
%% subplot(2,4,8);
subplot(2,4,8);
plotBoxplot_subplot(8,actFcData,c3actExec,c4actExec,hcSubj,'LIPS Activation','');

keyboard
%%
% %%
% actVar = c4c1fc;
% 
% rfefAct=actFcData(:,actVar);
% pdData=rfefAct(pdSubj);
% pdAge=age(pdSubj);
% hcData=rfefAct(hcSubj);
% hcAge=age(hcSubj);
% outlierThresh=5*std(rfefAct);
% pdSubj(pdData>outlierThresh)=[];
% hcSubj(hcData>outlierThresh)=[];
% pdData=rfefAct(pdSubj);
% pdAge=age(pdSubj);
% hcData=rfefAct(hcSubj);
% hcAge=age(hcSubj);
% nrCond=8;
% for plotNr=1:nrCond
%    subplot(1,nrCond,plotNr);
%    pdRt=rtData(pdSubj,plotNr);
%    hcRt=rtData(hcSubj,plotNr);
%    %plot(pdRt,pdData,'b.',hcRt,hcData,'r.');
%    %xlim([500,1450]);
%    ylim([-1,3]);
%    hold on
%    pdMdl=fitlm(pdRt,pdData);
%    hPa=plotAdded(pdMdl);
%    delete(hPa(3));
%    hPa(2).Color='b';
%    hcMdl=fitlm(hcRt,hcData);
%    hPa=plotAdded(hcMdl);
%    delete(hPa(3));
%    hold off
%    legend('off');
%    [p,F]=compareSlopes([pdRt(:);hcRt(:)],[pdData(:);hcData(:)],[ones(numel(pdData),1);zeros(numel(hcData),1)]);
%    disp([corr(pdRt,pdData),corr(hcRt,hcData),p,F]);
% end
% % corr
% %   Exec
% subSel=1:46;[rho,pval]=corr(rtData(subSel,rtExec),actFcData(subSel,c2actExec))%r=0.4002,p=0.0058 *
% subSel=1:46;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c2actExec),age(subSel))%r=0.4501,p=0.0019 *
% subSel=1:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec))%r=-0.3399,p=0.0208
% %%   Incon
% subSel=1:46;[rho,pval]=corr(rtData(subSel,rtIncon),actFcData(subSel,c2actIncon))%r=0.4033,p=
% subSel=1:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actIncon))%r=-0.2290,p=
% % partialcorr
% subSel=1:46;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c2actExec),age(subSel))%r=0.4501,p=0.0019
% subSel=1:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.2819,p=0.0606
% 
% %% pd/hc
% % c2
% % Exec
% %   partialcorr
% subSel=1:25;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c2actExec),age(subSel))%r=0.4255,p=0.0382 *
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.1159,p=0.59
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),rtData(subSel,rtExec),age(subSel))%r=-0.02
% subSel=26:46;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c2actExec),age(subSel))%r=0.3057,p=0.19 *
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.1159,p=0.6
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),rtData(subSel,rtExec),age(subSel))%r=0.15
% %   corr
% subSel=1:25;[rho,pval]=corr(rtData(subSel,rtExec),actFcData(subSel,c2actExec))%r=0.4006,p=0.0472 *
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec))%r=-0.1507,p=0.47
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2c1fc),rtData(subSel,rtExec))%r=-0.03
% subSel=26:46;[rho,pval]=corr(rtData(subSel,rtExec),actFcData(subSel,c2actExec))%r=0.2783,p=0.22 *
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec))%r=-0.12,p=0.6
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),rtData(subSel,rtExec))%r=0.21
% %% Incon
% %   partialcorr
% subSel=1:25;[rho,pval]=partialcorr(rtData(subSel,rtIncon),actFcData(subSel,c2actIncon),age(subSel))%r=0.22,p=0.29
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actIncon),age(subSel))%r=-0.08
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),rtData(subSel,rtIncon),age(subSel))%r=-0.14
% subSel=26:46;[rho,pval]=partialcorr(rtData(subSel,rtIncon),actFcData(subSel,c2actIncon),age(subSel))%r=0.464,p=0.039
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actIncon),age(subSel))%r=0.05
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),rtData(subSel,rtIncon),age(subSel))%r=-0.03
% %   corr
% subSel=1:25;[rho,pval]=corr(rtData(subSel,rtIncon),actFcData(subSel,c2actIncon))%r=0.31,p=0.11
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actIncon))%r=-0.13
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2c1fc),rtData(subSel,rtIncon))%r=-0.17
% subSel=26:46;[rho,pval]=corr(rtData(subSel,rtIncon),actFcData(subSel,c2actIncon))%r=0.39,p=0.077
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actIncon))%r=0.01
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),rtData(subSel,rtIncon))%r=0.03
% 
% %% c4
% % Exec
% %   partialcorr
% subSel=1:25;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c4actExec),age(subSel))%r=0.17
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actExec),age(subSel))%r=-0.002 *
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),rtData(subSel,rtExec),age(subSel))%r=-0.1
% subSel=26:46;[rho,pval]=partialcorr(rtData(subSel,rtExec),actFcData(subSel,c4actExec),age(subSel))%r=0.1
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actExec),age(subSel))%r=-0.6862,p=0.0008 *
% figure,plot(actFcData(1:25,c4c1fc),actFcData(1:25,c4actExec),'b*',actFcData(26:46,c4c1fc),actFcData(26:46,c4actExec),'r*')
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),rtData(subSel,rtExec),age(subSel))%r=-0.04
% %   corr
% subSel=1:25;[rho,pval]=corr(rtData(subSel,rtExec),actFcData(subSel,c4actExec))%r=0.1
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actExec))%r=-0.03 **
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c4c1fc),rtData(subSel,rtExec))%r=-0.1
% subSel=26:46;[rho,pval]=corr(rtData(subSel,rtExec),actFcData(subSel,c4actExec))%r=0.1
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actExec))%r=-0.6483,p=0.0015 **
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c1fc),rtData(subSel,rtExec))%r=0.07
% %% Incon
% %   partialcorr
% subSel=1:25;[rho,pval]=partialcorr(rtData(subSel,rtIncon),actFcData(subSel,c4actIncon),age(subSel))%r=-0.06
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actIncon),age(subSel))%r=0.03
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),rtData(subSel,rtIncon),age(subSel))%r=-0.1
% subSel=26:46;[rho,pval]=partialcorr(rtData(subSel,rtIncon),actFcData(subSel,c4actIncon),age(subSel))%r=0.1
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actIncon),age(subSel))%r=-0.02
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),rtData(subSel,rtIncon),age(subSel))%r=-0.08
% %   corr
% subSel=1:25;[rho,pval]=corr(rtData(subSel,rtIncon),actFcData(subSel,c4actIncon))%r=0.02
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actIncon))%r=-0.06
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c4c1fc),rtData(subSel,rtIncon))%r=-0.24
% subSel=26:46;[rho,pval]=corr(rtData(subSel,rtIncon),actFcData(subSel,c4actIncon))%r=0.08
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c4actIncon))%r=-0.1
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c1fc),rtData(subSel,rtIncon))%r=0.02
% 
% %% c2-c4 activation
% subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2actExec),actFcData(subSel,c4actExec),age(subSel))%r=0.71,p=0.000098
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2actExec),actFcData(subSel,c4actExec),age(subSel))%r=0.28,p=0.22
% subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2actExec),actFcData(subSel,c4actExec))%r=0.6973,p=0.00010701
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2actExec),actFcData(subSel,c4actExec))%r=0.2868,p=0.2
% 
% %% c4c2 & c4c3
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c2fc),actFcData(subSel,c4actExec),age(subSel))%r=-0.23,p=0.3
% subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c3fc),actFcData(subSel,c4actExec),age(subSel))%r=-0.38,p=0.09
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c2fc),actFcData(subSel,c4actExec))%r=-0.2
% subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c3fc),actFcData(subSel,c4actExec))%r=-0.38,p=0.08
% 
% g1 = [true(25,1);false(21,1)];
% g2=~g1;
% X = [g2, age,actFcData(:,c4actExec).*g1, actFcData(:,c4actExec).*g2];
% s = regstats(actFcData(:,c4c1fc),X);
% [p,F] = linhyptest(s.beta, s.covb, 0, [0 0 0 1 -1], s.tstat.dfe)
% 
% c2actOnC4fc=actFcData(:,c2actExec)./(actFcData(:,c4c1fc).*actFcData(:,c1actExec));
% figure
% plot(c2actOnC4fc)
% 
% figure
% subplot(1,2,1);
% plot(actFcData(1:25,c4c1fc),actFcData(1:25,c4actExec),'b*',actFcData(26:46,c4c1fc),actFcData(26:46,c4actExec),'r*');
% subplot(1,2,2);
% plot(actFcData(1:25,c2actExec),actFcData(1:25,c4actExec),'b*',actFcData(26:46,c2actExec),actFcData(26:46,c4actExec),'r*');
% 
% figure
% subplot(2,2,1);
% hPa=plotAdded(fitlm(actFcData(1:25,c4c1fc),actFcData(1:25,c4actExec)));
% xlim([floor(min(actFcData(:,c4c1fc))),ceil(max(actFcData(:,c4c1fc)))]);
% ylim([floor(min(actFcData(:,c4actExec))),ceil(max(actFcData(:,c4actExec)))]);
% xlabel('');ylabel('PD');title('');
% [r,p]=corr(actFcData(1:25,c4c1fc),actFcData(1:25,c4actExec));
% legend({['r=',num2str(r,2),',p=',num2str(p,1)]},'box','off');
% hPa(1).Color='b';hPa(2).Color='b';hPa(3).Color='b';
% subplot(2,2,2);
% hPa=plotAdded(fitlm(actFcData(1:25,c2actExec),actFcData(1:25,c4actExec)));
% xlim([floor(min(actFcData(:,c2actExec))),ceil(max(actFcData(:,c2actExec)))]);
% ylim([floor(min(actFcData(:,c4actExec))),ceil(max(actFcData(:,c4actExec)))]);
% xlabel('');ylabel('');title('');
% [r,p]=corr(actFcData(1:25,c2actExec),actFcData(1:25,c4actExec));
% legend({['r=',num2str(r,2),',p=',num2str(p,1)]},'box','off','FontSize',11);
% hPa(1).Color='b';hPa(2).Color='b';hPa(3).Color='b';
% subplot(2,2,3);
% hPa=plotAdded(fitlm(actFcData(26:46,c4c1fc),actFcData(26:46,c4actExec)));
% xlim([floor(min(actFcData(:,c4c1fc))),ceil(max(actFcData(:,c4c1fc)))]);
% ylim([floor(min(actFcData(:,c4actExec))),ceil(max(actFcData(:,c4actExec)))]);
% xlabel('Striatal-RIPS connectivity');ylabel('HC');title('');
% [r,p]=corr(actFcData(26:46,c4c1fc),actFcData(26:46,c4actExec));
% legend({['r=',num2str(r,2),',p=',num2str(p,1)]},'box','off','FontSize',11);
% hPa(2).Color='r';hPa(1).MarkerEdgeColor='r';
% text(-15,2.5,'RIPS Activation','Rotation',90,'FontSize',14)
% subplot(2,2,4);
% hPa=plotAdded(fitlm(actFcData(26:46,c2actExec),actFcData(26:46,c4actExec)));
% xlim([floor(min(actFcData(:,c2actExec))),ceil(max(actFcData(:,c2actExec)))]);
% ylim([floor(min(actFcData(:,c4actExec))),ceil(max(actFcData(:,c4actExec)))]);
% xlabel('RFEF Activation');ylabel('');title('');
% [r,p]=corr(actFcData(26:46,c2actExec),actFcData(26:46,c4actExec));
% legend({['r=',num2str(r,2),',p=',num2str(p,1)]},'box','off');
% hPa(2).Color='r';hPa(1).MarkerEdgeColor='r';
% %%


%% c4c1fc/c2actExec
% Exec
%   partialcorr
subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.04
subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c4c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.01
%   corr
subSel=1:25;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c2actExec))%r=-0.1
subSel=26:46;[rho,pval]=corr(actFcData(subSel,c4c1fc),actFcData(subSel,c2actExec))%r=-0.02

% c2c1fc/c4actExec
% Exec
%   partialcorr
subSel=1:25;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c4actExec),age(subSel))%r=-0.009
subSel=26:46;[rho,pval]=partialcorr(actFcData(subSel,c2c1fc),actFcData(subSel,c2actExec),age(subSel))%r=-0.1
%   corr
subSel=1:25;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c4actExec))%r=-0.02
subSel=26:46;[rho,pval]=corr(actFcData(subSel,c2c1fc),actFcData(subSel,c4actExec))%r=-0.35,p=0.11
%%
%rt=rtData([pdSubj(:);hcSubj(:)],1:5);
rt=rtData([pdSubj(:);hcSubj(:)],rtIncon);
act=[pdData(:);hcData(:)];
g=[ones(numel(pdData),1);zeros(numel(hcData),1)];
[p,F]=compareSlopes(rt,act,g);
keyboard
% bar(5:5:20,mnEdata);
% hold on
% errorbar(x+2.5,mnEdata,seEdata,'ko');
% hold off
% 
% figure
% actFcData=csvread('actFcData.csv');
% edata=actFcData(:,[1,2,4,5,7,8,10,11]);
% edata={edata(1:26,:);edata(27:46,:)};
% boxplot(edata{1});
% ylim1=ylim;
% boxplot(edata{2});
% ylim2=ylim;
% ylims=[min([ylim1(1),ylim2(1)]),max([ylim1(2),ylim2(2)])];
% condGroup=repmat([zeros(46,1);ones(46,1)],4,1);
% siteGroup=reshape(repmat(1:4,2*46,1),[],1);
% subplot(1,11,1);
% boxplot(edata{1}(:,1:2),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis on; ylim(ylims);
% set(gca,'Color','none','XColor','none','TickLength',[0,0],'box','off','YColor','k');
% subplot(1,11,2);
% boxplot(edata{2}(:,1:2),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,4);
% boxplot(edata{1}(:,3:4),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,5);
% boxplot(edata{2}(:,3:4),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,7);
% boxplot(edata{1}(:,5:6),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,8);
% boxplot(edata{2}(:,5:6),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,10);
% boxplot(edata{1}(:,7:8),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,11,11);
% boxplot(edata{2}(:,7:8),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% hOutliers=findobj(gcf,'Tag','Outliers');
% set(hOutliers,'Visible','off');
% hWhiskers=findobj(gcf,'Tag','Whisker');
% YData=reshape([hWhiskers.YData],2,[]);
% ylims=[floor(min(YData(1,:))),ceil(max(YData(2,:)))];
% hAxes=findall(gcf,'type','axes');
% for hNr=1:numel(hAxes)
%     ylim(hAxes(hNr),ylims);
% end
% 
% figure
% actFcData=csvread('actFcData.csv');
% edata=actFcData(:,[1,2,4,5,7,8,10,11]);
% edata={edata(1:26,:);edata(27:46,:)};
% boxplot(edata{1});
% ylim1=ylim;
% boxplot(edata{2});
% ylim2=ylim;
% ylims=[min([ylim1(1),ylim2(1)]),max([ylim1(2),ylim2(2)])];
% condGroup=repmat([zeros(46,1);ones(46,1)],4,1);
% siteGroup=reshape(repmat(1:4,2*46,1),[],1);
% subplot(1,12,1);
% axis on; ylim(ylims);
% set(gca,'Color','none','XColor','none','TickLength',[0,0],'box','off');
% subplot(1,12,2);
% boxplot(edata{1}(:,1:2),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,3);
% boxplot(edata{2}(:,1:2),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,5);
% boxplot(edata{1}(:,3:4),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,6);
% boxplot(edata{2}(:,3:4),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,8);
% boxplot(edata{1}(:,5:6),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,9);
% boxplot(edata{2}(:,5:6),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,11);
% boxplot(edata{1}(:,7:8),'plotstyle','compact','notch','on','colors','bc','jitter',0);
% axis off; ylim(ylims);
% subplot(1,12,12);
% boxplot(edata{2}(:,7:8),'plotstyle','compact','notch','on','colors','rm','jitter',0);
% axis off; ylim(ylims);
% hOutliers=findobj(gcf,'Tag','Outliers');
% set(hOutliers,'Visible','off');
% hWhiskers=findobj(gcf,'Tag','Whisker');
% YData=reshape([hWhiskers.YData],2,[]);
% ylims=[floor(min(YData(1,:))),ceil(max(YData(2,:)))];
% hAxes=findall(gcf,'type','axes');
% for hNr=1:numel(hAxes)
%     ylim(hAxes(hNr),ylims);
% end
% 
