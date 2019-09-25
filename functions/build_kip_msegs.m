function [ Kip ] = build_kip_msegs( entry, exit, nhat, P_all, GP,nsegs,options )
%[ Kip ] = build_kip( entry, exit, nhat, P_all, GP )
%   Detailed explanation goes here

if nargin < 7
    options = 'on';
end
if strcmp(options,'on')
    disp('Building Kip')
    tic
end
[n_test, ~] = size(P_all);
[~,n_obs] = size(entry);

Kpi = zeros(n_test*3,n_obs);

sig_f = GP.sig_f;
M = 1./GP.M.^2;
a = GP.a;
b = GP.b;

i_res = 1000;

maxSet = 50;
R_set = 1:n_obs;

while(~isempty(R_set))
    inds = 1:min(maxSet,length(R_set));
    set = R_set(inds);
    n_s = length(set);
    R_set(inds) = [];
    
    Ltotal = zeros(1,length(set));
    for ss = 1:max(nsegs(set))
%     for ss = 1:1
        Isegs = (nsegs(set) >= ss);        % index of measurements that have at least this many segments
        segSet = set(Isegs);
        n_ss = length(segSet);
        % create a whole bunch of points that go from entry to exit
        Rx = linspace2D(entry(2*ss-1,segSet)',exit(2*ss-1,segSet)',i_res );
        Ry = linspace2D(entry(2*ss,segSet)',exit(2*ss,segSet)',i_res );
        
        d = exit(2*ss-1:2*ss,segSet)-entry(2*ss-1:2*ss,segSet);
        L = hypot(d(1,:),d(2,:));
        
        Ltotal(Isegs) = Ltotal(Isegs)+L;
        
        N1 = repmat(nhat(1,segSet).^2,n_test,1,i_res);
        N12 = repmat(2*nhat(1,segSet).*nhat(2,segSet),n_test,1,i_res);
        N2 = repmat(nhat(2,segSet).^2,n_test,1,i_res);
        
        A = bsxfun(@minus,P_all(:,1),Rx(:)');
        A = reshape(A,n_test,n_ss,[]);
        B = bsxfun(@minus,P_all(:,2),Ry(:)');
        B = reshape(B,n_test,n_ss,[]);
        
        Kg = sig_f^2*exp(-(A.^2*M(1)+B.^2*M(2))/2);
        
        
        Kgdx2dx2  = M(1)^2*(A.^4*M(1)^2 - 6*M(1)*A.^2 + 3).*Kg;
        Kgdy2dy2 = M(2)^2*(B.^4*M(2)^2 - 6*M(2)*B.^2 + 3).*Kg;
        Kgdx2dy2 = M(1)*M(2)*(1-M(1)*A.^2).*(1-M(2)*B.^2).*Kg;
        Kgdxdy3 = M(1)*M(2)^2*A.*B.*(M(1)*B.^2-3).*Kg;
        Kgdx3dy = M(2)*M(1)^2*B.*A.*(M(2)*A.^2-3).*Kg;
        
        K11 = Kgdx2dx2 - 2*a/b*Kgdx2dy2 + (a/b)^2*Kgdy2dy2;
        K22 = (a+b)^2/b^2*Kgdx2dy2;
        K33 = (a/b)^2*Kgdx2dx2 - 2*a/b*Kgdx2dy2 + Kgdy2dy2;
        K12 = (a+b)/b*Kgdx3dy - a*(a+b)/b^2*Kgdxdy3;
        K13 = -a/b*Kgdx2dx2 + (a^2+b^2)/b^2*Kgdx2dy2 - a/b*Kgdy2dy2;
        K23 = (a+b)/b*Kgdxdy3 - a*(a+b)/b^2*Kgdx3dy;
        
        % Kpi1 = trapz(Kpi1,3).*(repmat(l',length(u_all),1)/(i_res-1)); %
        % calculates correct integral but then you need to divide by l for the
        % measurement so
        Kpi1= N1.*K11 + N12.*K12 + N2.*K13;
        Kpi1 = trapz(Kpi1,3)./(i_res-1);
        Kpi2 = N1.*K12 + N12.*K22 + N2.*K23;
        Kpi2 = trapz(Kpi2,3)./(i_res-1);
        Kpi3 = N1.*K13 + N12.*K23 + N2.*K33;
        Kpi3 = trapz(Kpi3,3)./(i_res-1);
        
        Kpi(1:3:end,segSet) = Kpi(1:3:end,segSet)+ Kpi1.*L;
        Kpi(2:3:end,segSet) = Kpi(2:3:end,segSet)+ Kpi2.*L;
        Kpi(3:3:end,segSet) = Kpi(3:3:end,segSet)+ Kpi3.*L;
    end
    Kpi(1:3:end,segSet) = Kpi(1:3:end,segSet)./Ltotal;
    Kpi(2:3:end,segSet) = Kpi(2:3:end,segSet)./Ltotal;
    Kpi(3:3:end,segSet) = Kpi(3:3:end,segSet)./Ltotal;
    if strcmp(options,'on')
        toc
    end
end
Kip = Kpi';

end

