function x = lsi(E,f,G,h);
% x = lsi(E,f,G,h)
% Least squares with inequality constraints.
% Solves the following problem:
%       min || Ex - f || subject to Gx >= h
% Reference: Lawson and Hanson (1974).
[Q,R] = qr(E);
[m,n] = size(R);
   Q1 = Q(:,1:n);
    R = R(1:n,:);
   RR = inv(R);
   GG = G * RR;
   f1 = Q1' * f;
   hh = h - GG * f1;
    z = ldp(GG,hh);
    x = RR*(z+f1);
