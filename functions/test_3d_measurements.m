%% test 1
Fxx = @(x,y,z) 1;
Fyy = @(x,y,z) 0;
Fzz = @(x,y,z) 0;
Fxy = @(x,y,z) 0;
Fxz = @(x,y,z) 0;
Fyz = @(x,y,z) 0;

nhat = [1;0;0];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y==1);

%% test 2
Fxx = @(x,y,z) 0;
Fyy = @(x,y,z) 2;
Fzz = @(x,y,z) 0;
Fxy = @(x,y,z) 0;
Fxz = @(x,y,z) 0;
Fyz = @(x,y,z) 0;

nhat = [0;1;0];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y==2);
%% test 3

Fxx = @(x,y,z) 0;
Fyy = @(x,y,z) 0;
Fzz = @(x,y,z) 3;
Fxy = @(x,y,z) 0;
Fxz = @(x,y,z) 0;
Fyz = @(x,y,z) 0;

nhat = [0;0;1];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y==3);

%% test 4

Fxx = @(x,y,z) 0;
Fyy = @(x,y,z) 0;
Fzz = @(x,y,z) 0;
Fxy = @(x,y,z) 1;
Fxz = @(x,y,z) 0;
Fyz = @(x,y,z) 0;

nhat = [1/sqrt(2);1/sqrt(2);0];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y<1+eps && y>1-eps);

%% test 5

Fxx = @(x,y,z) 0;
Fyy = @(x,y,z) 0;
Fzz = @(x,y,z) 0;
Fxy = @(x,y,z) 0;
Fxz = @(x,y,z) 0;
Fyz = @(x,y,z) 1;

nhat = [0;1/sqrt(2);1/sqrt(2)];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y<1+eps && y>1-eps);

%% test 6

Fxx = @(x,y,z) 0;
Fyy = @(x,y,z) 0;
Fzz = @(x,y,z) 0;
Fxy = @(x,y,z) 0;
Fxz = @(x,y,z) 1;
Fyz = @(x,y,z) 0;

nhat = [1/sqrt(2);0;1/sqrt(2)];
entry = [0;0;0];
exit = entry + nhat;

y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y<1+eps && y>1-eps);

%% test 7

t = 6e-3;               % 6mm thickness             
h = 10e-3;              % 10mm height
l = 20e-3;              % 20mm length
E = 200e9;              % 200GPa young's modulus
nu = 0.3;                % poissons ratio
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

nhat = [1;0;0];
entry = [0;h/4;0];
exit = entry + nhat*l;

y_check = make_measurements( entry, exit, nhat, @(x,y)Fxx(x,y,0), @(x,y)Fxy(x,y,0),@(x,y)Fyy(x,y,0));
y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y<y_check+eps && y>y_check-eps);

%% test 8

t = 6e-3;               % 6mm thickness             
h = 10e-3;              % 10mm height
l = 20e-3;              % 20mm length
E = 200e9;              % 200GPa young's modulus
nu = 0.3;                % poissons ratio
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

nhat = [1/sqrt(2);1/sqrt(2);0];
entry = [0;0;0];
exit = entry + nhat*l/4;

y_check = make_measurements( entry, exit, nhat, @(x,y)Fxx(x,y,0), @(x,y)Fxy(x,y,0),@(x,y)Fyy(x,y,0));
y = make_measurements_3D( entry, exit, nhat, Fxx, Fxy, Fxz,Fyy,Fyz,Fzz);

assert(y<y_check+eps && y>y_check-eps);
