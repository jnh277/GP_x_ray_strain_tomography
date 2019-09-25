function [ y ] = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz)
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
N11 = repmat(nhat(1,:)'.^2,1,i_res);
N12 = repmat(2*nhat(1,:)'.*nhat(2,:)',1,i_res);
N13 = repmat(2*nhat(1,:)'.*nhat(3,:)',1,i_res);
N22 = repmat(nhat(2,:)'.^2,1,i_res);
N23 = repmat(2*nhat(2,:)'.*nhat(3,:)',1,i_res);
N33 = repmat(nhat(3,:)'.^2,1,i_res);
% do math
V =  N11.*Vxx+N12.*Vxy+N13.*Vxz+N22.*Vyy+N23.*Vyz+N33.*Vzz;
% S = sqrt((Rx-Rx(:,1)).^2 + (Ry-Ry(:,1)).^2);
% L = S(:,end) - S(:,1);
% ds = L/i_res;
% integrate
y = trapz(1:i_res,V,2)/(i_res-1);

end

