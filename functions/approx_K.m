function [ Phi, Phi_pred, invLambda, sq_lambda, diagLambda ] = approx_K( nhat,entry, exit,P_all,m_basis, GP )
%[ Phi, Phi_pred, invLambda ] = approx_K( nhat,entry, exit,m_basis, GP )
%   Detailed explanation goes here

m_1 = m_basis; m_2 = m_basis;
sig_f = GP.sig_f;
M = GP.M;

a = GP.a;
b = GP.b;
A=a/b; B=-(a+b)/b;              % implementational form
L = sqrt(sum((exit-entry).^2))';

[n_test,~] = size(P_all);

mm1 = repmat(1:m_1,m_2,1);      % these basis functions are placed in a rectangle area which could be improved
mm2 = repmat([1:m_2]',1,m_1);
insideEllipse = mm1.^2/m_1^2 + mm2.^2/m_2^2 < 1;
mm1 = mm1(insideEllipse);
mm2 = mm2(insideEllipse);
MM = [mm1(:)'; mm2(:)'];   % all the basis functions, basis functions to change across columns
[~,mm_adj] = size(MM);

% extract info about lines
nx = nhat(1,:)';        % nx of every integral
ny = nhat(2,:)';        % ny of every integral
x0 = entry(1,:)';        % entry location of every integral
y0 = entry(2,:)';        % exit location of every integral


dlambda_x = 3.5/M(1)/m_1;
dlambda_y = 3.5/M(2)/m_2;

lambda_x = MM(1,:)*dlambda_x;
lambda_y = MM(2,:)*dlambda_y;
Lx = pi/2/dlambda_x;
Ly = pi/2/dlambda_y;

% lambda_x = MM(1,:)*pi/Lx/2;
% lambda_y = MM(2,:)*pi/Ly/2;
lambdaX = bsxfun(@times,nx,lambda_x);
lambdaY = bsxfun(@times,ny,lambda_y);
lambda_plus = lambdaX+lambdaY;
lambda_minus = lambdaX-lambdaY;
Bx = bsxfun(@times,x0+Lx,lambda_x);     
By = bsxfun(@times,y0+Ly,lambda_y);
B_plus = Bx+By;
B_minus = Bx-By;

s0 = zeros(size(L));
sf = L;

% BUILD PHI
Cs = bsxfun(@times,nx.^2,lambda_x.^2-A*lambda_y.^2) + bsxfun(@times,ny.^2,lambda_y.^2-A*lambda_x.^2);       
C0 = 2*B*bsxfun(@times,nx.*ny,lambda_x.*lambda_y);
I = Cs.*( (sin(lambda_minus.*repmat(sf,1,mm_adj)+B_minus)-sin(lambda_minus.*repmat(s0,1,mm_adj)+B_minus))./lambda_minus +...
                            (-sin(lambda_plus.*repmat(sf,1,mm_adj)+B_plus)+sin(lambda_plus.*repmat(s0,1,mm_adj)+B_plus))./lambda_plus )+...
                            C0.*( (sin(lambda_minus.*repmat(sf,1,mm_adj)+B_minus)-sin(lambda_minus.*repmat(s0,1,mm_adj)+B_minus))./lambda_minus +...
                            (sin(lambda_plus.*repmat(sf,1,mm_adj)+B_plus)-sin(lambda_plus.*repmat(s0,1,mm_adj)+B_plus))./lambda_plus );

Phi = I./repmat(L,1,mm_adj)/2/sqrt(Lx*Ly);

% BUILD PHI_pred
dx2 = -lambda_x.^2.*(1/sqrt(Ly*Lx)).*sin(bsxfun(@times,P_all(:,1)+Lx,lambda_x)).*sin(bsxfun(@times,P_all(:,2)+Ly,lambda_y));
dy2 = -lambda_y.^2.*(1/sqrt(Ly*Lx)).*sin(bsxfun(@times,P_all(:,1)+Lx,lambda_x)).*sin(bsxfun(@times,P_all(:,2)+Ly,lambda_y));
dxdy = lambda_x.*lambda_y.*(1/sqrt(Ly*Lx)).*cos(bsxfun(@times,P_all(:,1)+Lx,lambda_x)).*cos(bsxfun(@times,P_all(:,2)+Ly,lambda_y));


Phi_pred = NaN(3*n_test,mm_adj);
Phi_pred(1:3:end,:) = A*dy2-dx2;
Phi_pred(2:3:end,:) = B*dxdy;
Phi_pred(3:3:end,:) = A*dx2-dy2;


% spectral density
sq_lambda = [lambda_x'   lambda_y']; % angular frequencies
invLambda = 1./(2*pi*M(1)*M(2)*exp(-0.5*(M(1)^2*sq_lambda(:,1).^2+... 
    M(2)^2*sq_lambda(:,2).^2))*sig_f^2);
diagLambda = (2*pi*M(1)*M(2)*exp(-0.5*(M(1)^2*sq_lambda(:,1).^2+... 
    M(2)^2*sq_lambda(:,2).^2))*sig_f^2);

end

