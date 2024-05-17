function [U,V]=ed(S);
% [u,v]=ed(s);
%
% Gives SORTED eigenvalues in v, associated NORMALIZED eigenvectors in u
%
% NOTE: Giving only "ed(s)" produces EIGENVECTORS only 
% To get ONLY sorted EIGENVALUES use "eiv"
%
% The rank is determined numerically, dropping eigenvalues
% that are less then a predetermined tolerance (1e-9).
% Associated eigenvectors are eliminated.

tol=eps;
[U,V]=eig(S);
V=diag(V);
i=abs(V)>tol;
U=U(:,i);
V=V(i);
p=length(V);
[V,i]=sort(V);
j=p:-1:1;
V=V(j);
V=diag(V);
U=U(:,i(j));
[n,m]=size(U);
d=(sum(U.^2)).^.5;
U=U./(ones(n,1)*d);