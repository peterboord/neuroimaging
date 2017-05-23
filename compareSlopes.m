function [p,F]=compareSlopes(x,y,g1)

g2=~g1;
nrCols=size(x,2);
disp(size(x,1));
X=g2;
H=[];
for colNr=1:nrCols
    X = [X, x(:,colNr).*g1, x(:,colNr).*g2]; %#ok<AGROW>
    if isempty(H)
        H=[0 0 1 -1];
    else
        H=[[H,zeros(size(H,1),2)];[zeros(1,size(H,2)),1,-1]]; 
    end
end
s = regstats(y,X);
[p,F] = linhyptest(s.beta, s.covb, zeros(size(H,1),1), H, s.tstat.dfe);
