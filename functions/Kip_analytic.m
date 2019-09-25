function [Kip] = Kip_analytic( entry, exit, nhat, P_all, GP )
%[Kip] = Kip_analytic( entry, exit, nhat, P_all, GP )
%   Detailed explanation goes here


l1 = GP.M(1);l2 = GP.M(2);
sig_f = GP.sig_f;
a = GP.a;
b = GP.b;

[n_p,~] = size(P_all); 
[~,n_l] = size(entry);

L = hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));

x1 = L';

% converting all the points into line coordinates
px = bsxfun(@minus,P_all(:,1)',entry(1,:)');
py = bsxfun(@minus,P_all(:,2)',entry(2,:)');

Px_l = repmat(nhat(1,:)',1,n_p).*px + repmat(nhat(2,:)',1,n_p).*py;
Py_l = -repmat(nhat(2,:)',1,n_p).*px + repmat(nhat(1,:)',1,n_p).*py;


% for the entry point, ie start of integral
A = -Px_l;
B = -Py_l;

K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);

d3dx3 = (3/l1^4*A - 1/l1^6*A.^3).*K;
d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2).*K;
        
intd4dy4 = (B.^4/l2^8-6/l2^6*B.^2+3/l2^4).*exp(-B.^2/2/l2^2).*erf(A/sqrt(2)/l1).*l1.*sqrt(pi/2).*sig_f^2;
d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B).*K;
d3dy3 = (3/l2^4*B - 1/l2^6*B.^3).*K;
        
G1s  = d3dx3 -2*a/b*d3dxdy2 + (a/b)^2*intd4dy4;
G2s = (a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3;
G3s = -(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2 -(a/b)*intd4dy4;

% for the end of integral
A = repmat(x1,1,n_p)-Px_l;

K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);
d3dx3 = (3/l1^4*A - 1/l1^6*A.^3).*K;
d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2).*K;
% intd4dy4 = -sqrt(pi/2)*sig_f^2*l1*(1/l2^8*B.^4 - 6/l2^6*B.^2 + 3/l2^4).*...
%             exp(-B.^2/2/l2^2).*erf(-A/sqrt(2)/l1);
intd4dy4 = (B.^4/l2^8-6/l2^6*B.^2+3/l2^4).*exp(-B.^2/2/l2^2).*erf(A/sqrt(2)/l1).*l1.*sqrt(pi/2).*sig_f^2;
d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B).*K;
d3dy3 = (3/l2^4*B - 1/l2^6*B.^3).*K;

G1f  = d3dx3 -2*a/b*d3dxdy2 + (a/b)^2*intd4dy4;
G2f = (a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3;
G3f = -(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2 -(a/b)*intd4dy4;


% elements of Kip in line coordinates
Kip1 = (G1f-G1s)./repmat(x1,1,n_p);
Kip2 = (G2f-G2s)./repmat(x1,1,n_p);
Kip3 = (G3f-G3s)./repmat(x1,1,n_p);

% Kip_lc = [Kip1 Kip2 Kip3];

% Transformation matri to cartesian coordiantes
n1 = nhat(1,:)';
n2 = nhat(2,:)';
T11 = repmat(n1.^2,1,n_p);
T12 = repmat(-2*n1.*n2,1,n_p);
T13 = repmat(n2.^2,1,n_p);
T21 = repmat(n2.*n1,1,n_p);
T22 = repmat(n1.^2 - n2.^2,1,n_p);
T23 = repmat(-n1.*n2,1,n_p);
T31 = repmat(n2.^2,1,n_p);
T32 = repmat(2*n1.*n2,1,n_p);
T33 = repmat(n1.^2,1,n_p);


Kip = NaN(n_l,n_p*3);
Kip(:,1:3:end) = [Kip1.*T11+Kip2.*T12+Kip3.*T13];
Kip(:,2:3:end) = [Kip1.*T21+Kip2.*T22+Kip3.*T23];
Kip(:,3:3:end) = [Kip1.*T31+Kip2.*T32+Kip3.*T33];


end

