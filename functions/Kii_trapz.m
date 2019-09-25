function [Kii] = Kii_trapz( entry, exit,nhat, n_obs, GP )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

l1 = GP.M(1);l2 = GP.M(2);
sig_f = GP.sig_f;
a = GP.a;
b = GP.b;

% trapz parameters
i_res = 10000;
R_set = 1:n_obs;
maxSet = 100;

% stuff for line s
L = hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));
x1 = L';


Kii = NaN(n_obs,n_obs);

while(~isempty(R_set))
    inds = 1:min(maxSet,length(R_set));
    set = R_set(inds);
    n_s = length(set);
    R_set(inds) = [];
    
    % points along the line in cartesian
    Rx = linspace2D(entry(1,set)',exit(1,set)',i_res );
    Ry = linspace2D(entry(2,set)',exit(2,set)',i_res );
    
    px = reshape(bsxfun(@minus,Rx(:)',entry(1,:)'),n_obs,n_s,i_res);
    py = reshape(bsxfun(@minus,Ry(:)',entry(2,:)'),n_obs,n_s,i_res);
    
    Px_l = repmat(nhat(1,:)',1,n_s,i_res).*px + repmat(nhat(2,:)',1,n_s,i_res).*py;
    Py_l = -repmat(nhat(2,:)',1,n_s,i_res).*px + repmat(nhat(1,:)',1,n_s,i_res).*py;
    
    
    n1_l = bsxfun(@times,nhat(1,:)',nhat(1,set))+bsxfun(@times,nhat(2,:)',nhat(2,set));
    n2_l = -bsxfun(@times,nhat(2,:)',nhat(1,set)) + bsxfun(@times,nhat(1,:)',nhat(2,set));
    
    N1_l = n1_l.^2;
    N12_l = 2*n1_l.*n2_l;
    N2_l = n2_l.^2;


    %% for the beginning of line1
    A = -Px_l;
    B = -Py_l;
    
    K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);
    
    d3dx3 = (3/l1^4*A - 1/l1^6*A.^3);
    d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2);
    intd4dy4 = -sqrt(pi/2)*sig_f^2*l1*(1/l2^8*B.^4 - 6/l2^6*B.^2 + 3/l2^4).*...
        exp(-B.^2/2/l2^2).*erf(-A/sqrt(2)/l1);
    d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B);
    d3dy3 = (3/l2^4*B - 1/l2^6*B.^3);
    
    G1s  = (d3dx3 -2*a/b*d3dxdy2).*K + (a/b)^2*intd4dy4;
    G2s = ((a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3).*K;
    G3s = (-(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2).*K -(a/b)*intd4dy4;

    %% for the end of integral
    A = repmat(x1,1,n_s,i_res)-Px_l;
    
    K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);
    d3dx3 = (3/l1^4*A - 1/l1^6*A.^3);
    d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2);
    intd4dy4 = -sqrt(pi/2)*sig_f^2*l1*(1/l2^8*B.^4 - 6/l2^6*B.^2 + 3/l2^4).*...
        exp(-B.^2/2/l2^2).*erf(-A/sqrt(2)/l1);
    d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B);
    d3dy3 = (3/l2^4*B - 1/l2^6*B.^3);
    
    G1f  = (d3dx3 -2*a/b*d3dxdy2).*K + (a/b)^2*intd4dy4;
    G2f = ((a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3).*K;
    G3f = (-(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2).*K -(a/b)*intd4dy4;
    
    % elements of Kip in line 1s coordinates
    Kip1 = repmat(1./x1,1,n_s,i_res).*(G1f-G1s);
    Kip2 = repmat(1./x1,1,n_s,i_res).*(G2f-G2s);
    Kip3 = repmat(1./x1,1,n_s,i_res).*(G3f-G3s);
    
    
    Kii_T = N1_l.*Kip1 + N12_l.*Kip2 + N2_l.*Kip3;
    Kii(:,set) = trapz(Kii_T,3)/(i_res-1);
    
    disp(['Remaining = ' num2str(length(R_set))])
end

% Kii = (Kii+Kii')/2;
end
    
