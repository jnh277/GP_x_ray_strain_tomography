function [V,F] = refine_mesh(V,F,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:n 
    [nv,~] = size(V);
    [nf,~] = size(F);

    cx = mean(reshape(V(F',1),3,[]));
    cy = mean(reshape(V(F',2),3,[]));
    cz = mean(reshape(V(F',3),3,[]));

    V = [V;cx' cy' cz'];
    F = [F(:,1), F(:,2), (nv+1:nv+nf)';
        F(:,1), (nv+1:nv+nf)', F(:,3);
        (nv+1:nv+nf)', F(:,2), F(:,3)];
end

end

