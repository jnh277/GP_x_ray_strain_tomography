function [ Kpp ] = build_Kpp( P_all, GP,options )
%[ Kpp ] = build_Kpp( P_all, GP )
%   Detailed explanation goes here
if nargin < 3
    options = true;
end



sig_f = GP.sig_f;
M = 1./GP.M.^2;
a = GP.a;
b = GP.b;



[n_test, ~] = size(P_all);
if options
    disp('Building Kpp')
    tic
end

A = bsxfun(@minus,P_all(:,1),P_all(:,1)');
B = bsxfun(@minus,P_all(:,2),P_all(:,2)');

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

Kpp = zeros(3*n_test,3*n_test);

Kpp(1:3:end,1:3:end) = K11;
Kpp(1:3:end,2:3:end) = K12;
Kpp(1:3:end,3:3:end) = K13;
Kpp(2:3:end,1:3:end) = K12;
Kpp(2:3:end,2:3:end) = K22;
Kpp(2:3:end,3:3:end) = K23;
Kpp(3:3:end,1:3:end) = K13;
Kpp(3:3:end,2:3:end) = K23;
Kpp(3:3:end,3:3:end) = K33;
if options
    toc
end

end

