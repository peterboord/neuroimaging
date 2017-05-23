function [y,w] = pca_ica(x)

n=size(x,1);
[E,D]=eig(cov(x'));
v=E*D^(-0.5)*E'*x;
z=repmat(sqrt(sum(v.^2)),n,1).*v;
[EE,DD]=eig(cov(z'));
y=EE'*v;
w=EE';
end