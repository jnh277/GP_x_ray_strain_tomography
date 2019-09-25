function [ yhist,hist_strains ] = make_measurements_histogram( entry, exit, nhat, Fxx, Fxy, Fyy)
%[ y ] = make_measurements( entry, exit, nhat, Fxx, Fxy, Fyy)
%   Detailed explanation goes here
i_res = 100000;
% create a whole bunch of points that go from entry to exit
Rx = linspace2D(entry(1,:)',exit(1,:)',i_res );
Ry = linspace2D(entry(2,:)',exit(2,:)',i_res );
% evaluate the strain fields at these points
Vxx = Fxx(Rx,Ry);
Vxy = Fxy(Rx,Ry);
Vyy = Fyy(Rx,Ry);
N1 = repmat(nhat(1,:)'.^2,1,i_res);
N12 = repmat(2*nhat(1,:)'.*nhat(2,:)',1,i_res);
N2 = repmat(nhat(2,:)'.^2,1,i_res);
% do math
V =  N1.*Vxx+N12.*Vxy+N2.*Vyy;

centers = linspace(min(V(:)),max(V(:)),1000)

counts = hist(V',centers)


S = sqrt((Rx-Rx(:,1)).^2 + (Ry-Ry(:,1)).^2);
L = S(:,end) - S(:,1);
% ds = L/i_res;
% integrate
yhist = counts/1000./repmat(L',1000,1);
hist_strains = centers;

end

