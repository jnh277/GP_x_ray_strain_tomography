function [ y ] = make_measurements_xray_3D( entry, exit, kappa, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz)
%y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz)
%   Detailed explanation goes here
i_res = 200;
% create a whole bunch of points that go from entry to exit
Rx = linspace2D(entry(1,:)',exit(1,:)',i_res );
Ry = linspace2D(entry(2,:)',exit(2,:)',i_res );
Rz = linspace2D(entry(3,:)',exit(3,:)',i_res );
% evaluate the strain fields at these points
Vxx = Fxx(Rx,Ry,Rz);
Vxy = Fxy(Rx,Ry,Rz);
Vxz = Fxz(Rx,Ry,Rz);
Vyy = Fyy(Rx,Ry,Rz);
Vyz = Fyz(Rx,Ry,Rz);
Vzz = Fzz(Rx,Ry,Rz);
K11 = repmat(kappa(1,:)'.^2,1,i_res);
K12 = repmat(2*kappa(1,:)'.*kappa(2,:)',1,i_res);
K13 = repmat(2*kappa(1,:)'.*kappa(3,:)',1,i_res);
K22 = repmat(kappa(2,:)'.^2,1,i_res);
K23 = repmat(2*kappa(2,:)'.*kappa(3,:)',1,i_res);
K33 = repmat(kappa(3,:)'.^2,1,i_res);
% do math
V =  K11.*Vxx+K12.*Vxy+K13.*Vxz+K22.*Vyy+K23.*Vyz+K33.*Vzz;
% S = sqrt((Rx-Rx(:,1)).^2 + (Ry-Ry(:,1)).^2);
% L = S(:,end) - S(:,1);
% ds = L/i_res;
% integrate
y = trapz(1:i_res,V,2)/(i_res-1);

end

