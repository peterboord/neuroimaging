function [coeff,score,explained,mu]=fastEig(X,nrPc)
if nargin < 2
    nrPc=1;
end
mu = mean(X,2);
X = bsxfun(@minus,X,mu);
C = X*X';%X'*X;
[U,D] = eigs(C,nrPc);%[V,D] = eigs(C,nrPc);
diagD=diag(D);
clear C;
%eigVals(neighNrInd(vpNr),:)=diagD;
explained=100*diagD/sum(X(:).^2);
%U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
s = sqrt(abs(diagD));
coeff=U;
score=U'*X;
score=bsxfun(@(x,c)x./c,score,s');
end