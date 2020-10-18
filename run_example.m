clear all
clc

addpath ./functions
addpath ./intersection3D
addpath ./Beltrami_stress
addpath ./x_ray_strain_tomo

%% set number of projections here
num_angles = 3; % min 3 for it to work, too large a number will cause example to run very slowly

% can also adjust the number of segments from teh debye-scherrer ring to
% make measurements from by changing the below number
m_per_ray = 36;

% Not the added noise is random and so the results will vary on each run
% particularly with a low number of projections (measurements)

% enable or disable hyperparameter optimisation
HP_optimisation = false;        % change this to true to run, otherwise will use prefound values
%% define parameters of cantilever beam
t = 6e-3;               % 6mm thickness             
h = 10e-3;              % 10mm height
l = 20e-3;              % 20mm length
E = 200e9;              % 200GPa young's modulus
nu = 0.28;                % poissons ratio
Py = 2e3;               % 2KN load          in the negative y direction
Iyy = t*h^3/12;           % second moment of inertia
Pz = 1e3;               % 1KN load in the negative z direction
Izz = h*t^3/12;         % second moment of inertia

Fxx = @(x,y,z) Py.*(l-x).*y/E/Iyy + Pz.*(l-x).*z/E/Izz;
Fyy = @(x,y,z) -nu*Py.*(l-x).*y/E/Iyy + -nu*Pz.*(l-x).*z/E/Izz;
Fzz = @(x,y,z) -nu*Py.*(l-x).*y/E/Iyy + -nu*Pz.*(l-x).*z/E/Izz;
Fxy = @(x,y,z) -(1+nu)/E*Py*(h^2/4-y.^2)/2/Iyy;
Fxz = @(x,y,z) -(1+nu)/E*Pz*(t^2/4-z.^2)/2/Izz;
Fyz = @(x,y,z) 0.*x + 0.*y +0.*z;

clear E;
%% load geometry
CB = load('CB_3D.mat');

%% generate rays
disp('Simulating measurements.')
% so could probably use 0.1mm 
zangles = linspace(0,pi,num_angles+1);
zangles = zangles(1:end-1);
angles = [zeros(2,num_angles);zangles];     % only rotations about z are non zero
num_pixels = 40; 
width = 30e-3;

l = max(CB.V(:,1));

[ lines ] = gen_rays_3d(angles, num_pixels, width);
lines(1,:) = lines(1,:) + l/2;

%% find intersections
[entry_tmp,exit_tmp,nhat_tmp,L_tmp,yInds] = find_intersects_3D(CB,lines);

%% compute transverse measurement directions (have assumed that there are no lines in the dir


uhat = repmat([0;0;1],1,length(nhat_tmp));

angs = linspace(0,2*pi,m_per_ray+1);
angs = angs(1:end-1);
diff_theta = deg2rad(5);
nhat = repmat(nhat_tmp,1,m_per_ray);
entry = repmat(entry_tmp,1,m_per_ray);
exit = repmat(exit_tmp,1,m_per_ray);
L = repmat(L_tmp,1,m_per_ray);

% calucalte hte almost perpendicular measurement directions
kappa = nan(3,length(nhat_tmp)*m_per_ray);
for i = 1:m_per_ray
    for j = 1:length(nhat_tmp)
        % build a rotation about the line direction
        R = cos(angs(i))*eye(3)+sin(angs(i))*skew(nhat_tmp(:,j))+(1-cos(angs(i)))*(nhat_tmp(:,j)*nhat_tmp(:,j)');
        vhat = cross(uhat(:,j),nhat_tmp(:,j));
        R2 = cos(diff_theta)*eye(3)+sin(diff_theta)*skew(vhat)+(1-cos(diff_theta))*(vhat*vhat');
        
        kappa(:,(i-1)*length(nhat_tmp)+j) = R*uhat(:,j);
            
    end
end


%% simulate measurements
y = make_measurements_xray_3D( entry, exit, kappa, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

% add some noise
sig_m = 1e-4;
y = y + sig_m*randn(size(y));

%% GP stuff

num_basis = 10;
lx = 1;
ly = 1;
lz = 1;
covFunc = 'SE';

A.lx = lx;
A.ly = ly;
A.lz = lz;
A.sig_f = 01;
A.covFunc = covFunc;

B.lx = lx;
B.ly = ly;
B.lz = lz;
B.sig_f = 1;
B.covFunc = covFunc;

C.lx = lx;
C.ly = ly;
C.lz = lz;
C.sig_f = 1;
C.covFunc = covFunc;

D.lx = lx;
D.ly = ly;
D.lz = lz;
D.sig_f = 1;
D.covFunc = covFunc;

E.lx = lx;
E.ly = ly;
E.lz = lz;
E.sig_f = 1;
E.covFunc = covFunc;

F.lx = lx;
F.ly = ly;
F.lz = lz;
F.sig_f = 1;
F.covFunc = covFunc;
%% Hyperparameter optimisation
% Currently using
if HP_optimisation
    
%     tmp = repmat([10,1,1,1],6,1);tmp = tmp';
%     theta0 = tmp(:);
    % start point was found by running optimisation from several different
    % locations and picking the one that resulted in the lowest objective
    % evaluation
    % Using simulated annealing so as to not have to deal with multiple
    % start points to avoid local minima
    theta0 = [100;100;100;100;100;100;100;
        0.5;50;100;50;20;100;50;100;100;25;50;3;0.5;50;10;0.5;10];
    max_iters = 10;
    [LogL,thetaopt] = run_HP_optimisation(theta0,num_basis,A,B,C,D,E,F,nu,entry,exit,y,sig_m,kappa,max_iters);

else
    thetaopt = [157.095950425404;94.1423527889590;136.451203677699;126.616503019542;155.410751276087;148.380528993722;124.030121014359;0.582106994745147;58.0866137310982;92.5090638982722;68.5389582852935;23.5678621164255;95.6058730588661;79.7986594521732;87.1155121104460;95.3351837557044;34.5201905786993;40.4771011544049;3.66109303654358;0.296584941497878;60.3646650317464;10.9383516307489;0.282132376258170;8.18044187576777];
end

A.sig_f = thetaopt(1);
A.lx = thetaopt(2);
A.ly = thetaopt(3);
A.lz = thetaopt(4);

B.sig_f = thetaopt(5);
B.lx = thetaopt(6);
B.ly = thetaopt(7);
B.lz = thetaopt(8);

C.sig_f = thetaopt(9);
C.lx = thetaopt(10);
C.ly = thetaopt(11);
C.lz = thetaopt(12);

D.sig_f = thetaopt(13);
D.lx = thetaopt(14);
D.ly = thetaopt(15);
D.lz = thetaopt(16);

E.sig_f = thetaopt(17);
E.lx = thetaopt(18);
E.ly = thetaopt(19);
E.lz = thetaopt(20);

F.sig_f = thetaopt(21);
F.lx = thetaopt(22);
F.ly = thetaopt(23);
F.lz = thetaopt(24);
%% making Cz
disp('Performing reconstructions and plotting')
nx = 20; ny = 10;nz = 6;

xgv = linspace(0,l,nx);
ygv = linspace(-h/2,h/2,ny);
zgv = linspace(-t/2,t/2,nz);
[X,Y,Z] = meshgrid(xgv,ygv,zgv);

[Phi, SLambda,~,Phi_T] = beltrami_approx_xray(num_basis,A,B,C,D,E,F,nu,entry,exit,kappa,X(:),Y(:),Z(:));
[n,m]= size(Phi);
Gamma = [(ones(n,1)/sig_m).*Phi;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);
optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
vbasis = (linsolve(CZ,linsolve(CZ,Phi'*y/sig_m^2,optsT),opts));
%% Now reconstruct on the different cutting planes and plot
V = CB.V;

xp = l/5;
zp = 0;
xcut = V(:,1) > xp;
zcut = V(:,3) > zp;

V1 = V;
V1(xcut,1) = xp;

V2 = V;
V2(zcut,3) = zp;


az = 50.8000;
el = 33.2000;

tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;

figure(3)
clf
eff_ax(1)=axes;
p1 = patch('Faces',CB.F,'Vertices',V1,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
p2 = patch('Faces',CB.F,'Vertices',V2,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
camlight; camlight(-80,-10); 
view(eff_ax(1),[az el])
axis equal
title('Effective strain reconstruction error * 20','Interpreter','Latex','FontSize',20)
xlabel('x','Interpreter','Latex','FontSize',18)
ylabel('y','Interpreter','Latex','FontSize',18)
zlabel('z','Interpreter','Latex','FontSize',18)



figure(4)
clf
hyd_ax(1)=axes;
p1 = patch('Faces',CB.F,'Vertices',V1,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
p2 = patch('Faces',CB.F,'Vertices',V2,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
camlight; camlight(-80,-10); 
view(hyd_ax(1),[az el])
axis equal
title('Hydrostatic strain reconstruction error *20','Interpreter','Latex','FontSize',20)
xlabel('x','Interpreter','Latex','FontSize',18)
ylabel('y','Interpreter','Latex','FontSize',18)
zlabel('z','Interpreter','Latex','FontSize',18)



nx = 40; ny = nx/2;
xgv = linspace(xp,l,nx);
ygv = linspace(-h/2,h/2,ny);

[X,Y] = meshgrid(xgv,ygv);
Z = zp*ones(size(Y));

% construct kernel
[~, ~,~,Phi_T] = beltrami_approx_xray(num_basis,A,B,C,D,E,F,nu,entry,exit,kappa,X,Y,Z);
% f_approx = Phi_T*(linsolve(CZ,linsolve(CZ,Phi'*y/sig_m^2,optsT),opts));
f_approx = Phi_T*vbasis;

epsxx_pred = NaN(ny,nx);
epsyy_pred = NaN(ny,nx);
epszz_pred = NaN(ny,nx);
epsxy_pred = NaN(ny,nx);
epsxz_pred = NaN(ny,nx);
epsyz_pred = NaN(ny,nx);

epsxx_pred(:) = f_approx(1:6:end);
epsyy_pred(:) = f_approx(2:6:end);
epszz_pred(:) = f_approx(3:6:end);
epsxy_pred(:) = f_approx(4:6:end);
epsxz_pred(:) = f_approx(5:6:end);
epsyz_pred(:) = f_approx(6:6:end);



eps_hyd_pred = (epsxx_pred+epsyy_pred+epszz_pred)/3;
eps_eff_pred = sqrt(((epsxx_pred-eps_hyd_pred).^2 +(epsyy_pred-eps_hyd_pred).^2+...
                (epszz_pred-eps_hyd_pred).^2 + 2*(epsxy_pred.^2+epsxz_pred.^2 ...
                +epsyz_pred.^2))*2/3);
            
epsxx_true = NaN(ny,nx);
epsyy_true = NaN(ny,nx);
epszz_true = NaN(ny,nx);
epsxy_true = NaN(ny,nx);
epsxz_true = NaN(ny,nx);
epsyz_true = NaN(ny,nx);

epsxx_true(:) = Fxx(X(:),Y(:),Z(:));
epsyy_true(:) = Fyy(X(:),Y(:),Z(:));
epszz_true(:) = Fzz(X(:),Y(:),Z(:));
epsxy_true(:) = Fxy(X(:),Y(:),Z(:));
epsxz_true(:) = Fxz(X(:),Y(:),Z(:));
epsyz_true(:) = Fyz(X(:),Y(:),Z(:));


eps_hyd_true = (epsxx_true+epsyy_true+epszz_true)/3;
eps_eff_true = sqrt(((epsxx_true-eps_hyd_true).^2 +(epsyy_true-eps_hyd_true).^2+...
                (epszz_true-eps_hyd_true).^2 + 2*(epsxy_true.^2+epsxz_true.^2 ...
                +epsyz_true.^2))*2/3);

diff_eff = eps_eff_pred - eps_eff_true;
diff_hyd = eps_hyd_pred - eps_hyd_true;

            
min_eff = 1.85e-4;
max_eff = 0.0025;
min_hyd = -2.35e-4;
max_hyd = 4.3e-4;


            

figure(3)
eff_ax(2)=axes;
set(eff_ax(2),'color','none', ...
      'Position',get(eff_ax(1),'Position'), ...
      'Xlim',get(eff_ax(1),'Xlim'), ...
      'Ylim',get(eff_ax(1),'Ylim'), ...
      'Zlim',get(eff_ax(1),'Zlim'), ...
      'XTick',[], ...
      'YTick',[], ...
      'ZTick',[],...
      'Visible','off');
axis equal
view(eff_ax(2),[az el])
linkprop(eff_ax,{'view'});
% s = surface(X,Y,Z,eps_eff_pred);
error_scale = 20;
s = surface(X,Y,Z,diff_eff*error_scale);
set(s,'LineStyle','none')
% view([az el])

figure(4)
hyd_ax(2)=axes;
set(hyd_ax(2),'color','none', ...
      'Position',get(hyd_ax(1),'Position'), ...
      'Xlim',get(hyd_ax(1),'Xlim'), ...
      'Ylim',get(hyd_ax(1),'Ylim'), ...
      'Zlim',get(hyd_ax(1),'Zlim'), ...
      'XTick',[], ...
      'YTick',[], ...
      'ZTick',[],...
      'Visible','off');
axis equal
view(hyd_ax(2),[az el])
linkprop(eff_ax,{'view'});
% s = surface(X,Y,Z,eps_hyd_pred);
s = surface(X,Y,Z,diff_hyd*error_scale);
set(s,'LineStyle','none')

%% next cutting plane
nz = round(nx/3); ny = nx/2;
zgv = linspace(zp,t/2,nz);
ygv = linspace(-h/2,h/2,ny);

[Z,Y] = meshgrid(zgv,ygv);
X = xp*ones(size(Y));

[~, ~,~,Phi_T] = beltrami_approx_xray(num_basis,A,B,C,D,E,F,nu,entry,exit,kappa,X,Y,Z);
% f_approx = Phi_T*(linsolve(CZ,linsolve(CZ,Phi'*y/sig_m^2,optsT),opts));
f_approx = Phi_T*vbasis;

epsxx_pred = NaN(ny,nz);
epsyy_pred = NaN(ny,nz);
epszz_pred = NaN(ny,nz);
epsxy_pred = NaN(ny,nz);
epsxz_pred = NaN(ny,nz);
epsyz_pred = NaN(ny,nz);

epsxx_pred(:) = f_approx(1:6:end);
epsyy_pred(:) = f_approx(2:6:end);
epszz_pred(:) = f_approx(3:6:end);
epsxy_pred(:) = f_approx(4:6:end);
epsxz_pred(:) = f_approx(5:6:end);
epsyz_pred(:) = f_approx(6:6:end);

eps_hyd_pred = (epsxx_pred+epsyy_pred+epszz_pred)/3;
eps_eff_pred = sqrt(((epsxx_pred-eps_hyd_pred).^2 +(epsyy_pred-eps_hyd_pred).^2+...
                (epszz_pred-eps_hyd_pred).^2 + 2*(epsxy_pred.^2+epsxz_pred.^2 ...
                +epsyz_pred.^2))*2/3);
            
epsxx_true = NaN(ny,nz);
epsyy_true = NaN(ny,nz);
epszz_true = NaN(ny,nz);
epsxy_true = NaN(ny,nz);
epsxz_true = NaN(ny,nz);
epsyz_true = NaN(ny,nz);

epsxx_true(:) = Fxx(X(:),Y(:),Z(:));
epsyy_true(:) = Fyy(X(:),Y(:),Z(:));
epszz_true(:) = Fzz(X(:),Y(:),Z(:));
epsxy_true(:) = Fxy(X(:),Y(:),Z(:));
epsxz_true(:) = Fxz(X(:),Y(:),Z(:));
epsyz_true(:) = Fyz(X(:),Y(:),Z(:));


eps_hyd_true = (epsxx_true+epsyy_true+epszz_true)/3;
eps_eff_true = sqrt(((epsxx_true-eps_hyd_true).^2 +(epsyy_true-eps_hyd_true).^2+...
                (epszz_true-eps_hyd_true).^2 + 2*(epsxy_true.^2+epsxz_true.^2 ...
                +epsyz_true.^2))*2/3);

diff_eff = eps_eff_pred - eps_eff_true;
diff_hyd = eps_hyd_pred - eps_hyd_true;




figure(3)
s = surface(X,Y,Z,diff_eff*error_scale);
set(s,'LineStyle','none')
caxis([-max(abs([min_eff;max_eff])) max(abs([min_eff;max_eff]))])

colormap(cmap)
shading flat;

rectV = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1; 1 0 1]-0.5;
rectV(:,1) = rectV(:,1)*l+l/2;
rectV(:,2) = rectV(:,2)*h;
rectV(:,3) = rectV(:,3)*t;
rectF = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
Fshow = rectF([3 4 6],:);
p = patch('Vertices', rectV, 'Faces', Fshow);
set(p,'FaceAlpha',0,'LineStyle','-.')
line([xp xp],[-h/2 h/2],[0 0],'LineStyle','-.','Color','k')

figure(4)
s = surface(X,Y,Z,diff_hyd*error_scale);
set(s,'LineStyle','none')
caxis([-max(abs([min_hyd;max_hyd])) max(abs([min_hyd;max_hyd]))])
colormap(cmap)
shading flat;
line([xp xp],[-h/2 h/2],[0 0],'LineStyle','-.','Color','k')

rectV = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1; 1 0 1]-0.5;
rectV(:,1) = rectV(:,1)*l+l/2;
rectV(:,2) = rectV(:,2)*h;
rectV(:,3) = rectV(:,3)*t;
rectF = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
Fshow = rectF([3 4 6],:);

p = patch('Vertices', rectV, 'Faces', Fshow);
set(p,'FaceAlpha',0,'LineStyle','-.')

%% Plot true cutting planes
V = CB.V;

xp = l/5;
zp = 0;
xcut = V(:,1) > xp;
zcut = V(:,3) > zp;

V1 = V;
V1(xcut,1) = xp;

V2 = V;
V2(zcut,3) = zp;


az = 50.8000;
el = 33.2000;

tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;

figure(1)
clf
eff_ax(1)=axes;
p1 = patch('Faces',CB.F,'Vertices',V1,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
p2 = patch('Faces',CB.F,'Vertices',V2,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
camlight; camlight(-80,-10); 
view(eff_ax(1),[az el])
axis equal
title('True Effective Strain','Interpreter','Latex','FontSize',20)
xlabel('x','Interpreter','Latex','FontSize',18)
ylabel('y','Interpreter','Latex','FontSize',18)
zlabel('z','Interpreter','Latex','FontSize',18)



figure(2)
clf
hyd_ax(1)=axes;
p1 = patch('Faces',CB.F,'Vertices',V1,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
p2 = patch('Faces',CB.F,'Vertices',V2,'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
camlight; camlight(-80,-10); 
view(hyd_ax(1),[az el])
axis equal
title('True Hydrostatic Strain','Interpreter','Latex','FontSize',20)
xlabel('x','Interpreter','Latex','FontSize',18)
ylabel('y','Interpreter','Latex','FontSize',18)
zlabel('z','Interpreter','Latex','FontSize',18)



nx = 40; ny = nx/2;
xgv = linspace(xp,l,nx);
ygv = linspace(-h/2,h/2,ny);

[X,Y] = meshgrid(xgv,ygv);
Z = zp*ones(size(Y));



epsxx_pred = NaN(ny,nx);
epsyy_pred = NaN(ny,nx);
epszz_pred = NaN(ny,nx);
epsxy_pred = NaN(ny,nx);
epsxz_pred = NaN(ny,nx);
epsyz_pred = NaN(ny,nx);
% these are in fact the true strains but im too lazy to rename variables in the plotting code
epsxx_pred(:) = Fxx(X(:),Y(:),Z(:));  
epsyy_pred(:) = Fyy(X(:),Y(:),Z(:));
epszz_pred(:) = Fzz(X(:),Y(:),Z(:));
epsxy_pred(:) = Fxy(X(:),Y(:),Z(:));
epsxz_pred(:) = Fxz(X(:),Y(:),Z(:));
epsyz_pred(:) = Fyz(X(:),Y(:),Z(:));

eps_hyd_pred = (epsxx_pred+epsyy_pred+epszz_pred)/3;
eps_eff_pred = sqrt(((epsxx_pred-eps_hyd_pred).^2 +(epsyy_pred-eps_hyd_pred).^2+...
                (epszz_pred-eps_hyd_pred).^2 + 2*(epsxy_pred.^2+epsxz_pred.^2 ...
                +epsyz_pred.^2))*2/3);
            
        
            

figure(1)
eff_ax(2)=axes;
set(eff_ax(2),'color','none', ...
      'Position',get(eff_ax(1),'Position'), ...
      'Xlim',get(eff_ax(1),'Xlim'), ...
      'Ylim',get(eff_ax(1),'Ylim'), ...
      'Zlim',get(eff_ax(1),'Zlim'), ...
      'XTick',[], ...
      'YTick',[], ...
      'ZTick',[],...
      'Visible','off');
axis equal
view(eff_ax(2),[az el])
linkprop(eff_ax,{'view'});
s = surface(X,Y,Z,eps_eff_pred);
set(s,'LineStyle','none')


figure(2)
hyd_ax(2)=axes;
set(hyd_ax(2),'color','none', ...
      'Position',get(hyd_ax(1),'Position'), ...
      'Xlim',get(hyd_ax(1),'Xlim'), ...
      'Ylim',get(hyd_ax(1),'Ylim'), ...
      'Zlim',get(hyd_ax(1),'Zlim'), ...
      'XTick',[], ...
      'YTick',[], ...
      'ZTick',[],...
      'Visible','off');
axis equal
view(hyd_ax(2),[az el])
linkprop(eff_ax,{'view'});
s = surface(X,Y,Z,eps_hyd_pred);
set(s,'LineStyle','none')

% 
nz = round(nx/3); ny = nx/2;
zgv = linspace(zp,t/2,nz);
ygv = linspace(-h/2,h/2,ny);

[Z,Y] = meshgrid(zgv,ygv);
X = xp*ones(size(Y));


epsxx_pred = NaN(ny,nz); 
epsyy_pred = NaN(ny,nz);
epszz_pred = NaN(ny,nz);
epsxy_pred = NaN(ny,nz);
epsxz_pred = NaN(ny,nz);
epsyz_pred = NaN(ny,nz);

epsxx_pred(:) = Fxx(X(:),Y(:),Z(:));
epsyy_pred(:) = Fyy(X(:),Y(:),Z(:));
epszz_pred(:) = Fzz(X(:),Y(:),Z(:));
epsxy_pred(:) = Fxy(X(:),Y(:),Z(:));
epsxz_pred(:) = Fxz(X(:),Y(:),Z(:));
epsyz_pred(:) = Fyz(X(:),Y(:),Z(:));

eps_hyd_pred = (epsxx_pred+epsyy_pred+epszz_pred)/3;
eps_eff_pred = sqrt(((epsxx_pred-eps_hyd_pred).^2 +(epsyy_pred-eps_hyd_pred).^2+...
                (epszz_pred-eps_hyd_pred).^2 + 2*(epsxy_pred.^2+epsxz_pred.^2 ...
                +epsyz_pred.^2))*2/3);



figure(1)
s = surface(X,Y,Z,eps_eff_pred);
set(s,'LineStyle','none')
caxis([-max(abs([min_eff;max_eff])) max(abs([min_eff;max_eff]))])
colormap(cmap)
shading flat;

rectV = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1; 1 0 1]-0.5;
rectV(:,1) = rectV(:,1)*l+l/2;
rectV(:,2) = rectV(:,2)*h;
rectV(:,3) = rectV(:,3)*t;
rectF = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
Fshow = rectF([3 4 6],:);
p = patch('Vertices', rectV, 'Faces', Fshow);
set(p,'FaceAlpha',0,'LineStyle','-.')
line([xp xp],[-h/2 h/2],[0 0],'LineStyle','-.','Color','k')

figure(2)
s = surface(X,Y,Z,eps_hyd_pred);
set(s,'LineStyle','none')
caxis([-max(abs([min_hyd;max_hyd])) max(abs([min_hyd;max_hyd]))])
colormap(cmap)
shading flat;
line([xp xp],[-h/2 h/2],[0 0],'LineStyle','-.','Color','k')

rectV = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1; 1 0 1]-0.5;
rectV(:,1) = rectV(:,1)*l+l/2;
rectV(:,2) = rectV(:,2)*h;
rectV(:,3) = rectV(:,3)*t;
rectF = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
Fshow = rectF([3 4 6],:);

p = patch('Vertices', rectV, 'Faces', Fshow);
set(p,'FaceAlpha',0,'LineStyle','-.')

