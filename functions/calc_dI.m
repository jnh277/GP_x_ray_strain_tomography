function [dKii] = calc_dI(entry, exit, nhat, nhat_l2,entry_l2,s, GP )
    l1 = GP.M(1);l2 = GP.M(2);
    sig_f = GP.sig_f;
    a = GP.a;
    b = GP.b;
    
    n_p = length(s); 
    L = hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));
%     nhat = (exit-entry)./L;     % why is this error carrying through?????
    
    x1 = L';
    
    P_all = nhat_l2*s'+entry_l2;
    
    px = bsxfun(@minus,P_all(1,:),entry(1,:)');
	py = bsxfun(@minus,P_all(2,:),entry(2,:)');
    
    Px_l = repmat(nhat(1,:)',1,n_p).*px + repmat(nhat(2,:)',1,n_p).*py;
    Py_l = -repmat(nhat(2,:)',1,n_p).*px + repmat(nhat(1,:)',1,n_p).*py;
    
    n1_l = bsxfun(@times,nhat(1,:)',nhat_l2(1,:))+bsxfun(@times,nhat(2,:)',nhat_l2(2,:));
    n2_l = -bsxfun(@times,nhat(2,:)',nhat_l2(1,:)) + bsxfun(@times,nhat(1,:)',nhat_l2(2,:));
    
    N1_l = n1_l.^2;
    N12_l = 2*n1_l.*n2_l;
    N2_l = n2_l.^2;


    % for the entry point, ie start of integral
    A = -Px_l;
    B = -Py_l;

    K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);

    d3dx3 = (3/l1^4*A - 1/l1^6*A.^3).*K;
    d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2).*K;
    intd4dy4 = -sqrt(pi/2)*sig_f^2*l1*(1/l2^8*B.^4 - 6/l2^6*B.^2 + 3/l2^4).*...
                exp(-B.^2/2/l2^2).*erf(-A/sqrt(2)/l1);
%     intd4dy4 = (B.^4/l2^8-6/l2^6*B.^2+3/l2^4).*exp(-B.^2/2/l2^2).*erf(A/sqrt(2)/l1).*l1.*sqrt(pi/2);
    d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B).*K;
    d3dy3 = (3/l2^4*B - 1/l2^6*B.^3).*K;

    G1s  = d3dx3 -2*a/b*d3dxdy2 + (a/b)^2*intd4dy4;
    G2s = (a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3;
    G3s = -(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2 -(a/b)*intd4dy4;

    % for the end of integral
    A = repmat(x1,1,n_p)-Px_l;
    % B = -Py_l;

    K = sig_f^2*exp(-A.^2/2/l1^2-B.^2/2/l2^2);
    d3dx3 = (3/l1^4*A - 1/l1^6*A.^3).*K;
    d3dxdy2 = (1/l1^2/l2^2*A - 1/l1^2/l2^4*A.*B.^2).*K;
    intd4dy4 = -sqrt(pi/2)*sig_f^2*l1*(1/l2^8*B.^4 - 6/l2^6*B.^2 + 3/l2^4).*...
                exp(-B.^2/2/l2^2).*erf(-A/sqrt(2)/l1);
%     intd4dy4 = (B.^4/l2^8-6/l2^6*B.^2+3/l2^4).*exp(-B.^2/2/l2^2).*erf(A/sqrt(2)/l1).*l1.*sqrt(pi/2);
    d3dydx2 = (1/l1^2/l2^2*B - 1/l1^4/l2^2*A.^2.*B).*K;
    d3dy3 = (3/l2^4*B - 1/l2^6*B.^3).*K;

    G1f  = d3dx3 -2*a/b*d3dxdy2 + (a/b)^2*intd4dy4;
    G2f = (a+b)/b*d3dydx2 -a*(a+b)/b^2*d3dy3;
    G3f = -(a/b)*d3dx3 + (a^2+b^2)/b^2*d3dxdy2 -(a/b)*intd4dy4;

    % elements of Kip in line coordinates
    Kip1 = (G1f-G1s)./repmat(x1,1,n_p);
    Kip2 = (G2f-G2s)./repmat(x1,1,n_p);
    Kip3 = (G3f-G3s)./repmat(x1,1,n_p);
    
    dKii = N1_l.*Kip1 + N12_l.*Kip2 + N2_l.*Kip3;




end