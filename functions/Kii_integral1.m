function [Kii] = Kii_integral1( entry, exit,nhat, n_obs, GP,parallel )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

L = hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));
% nhat = (exit-entry)./L;
Kii = zeros(n_obs,n_obs);

parfor_progress(n_obs);
if ~parallel  
    for ii = 1:n_obs
        func = @(s) calc_dI(entry(:,ii:n_obs), exit(:,ii:n_obs), nhat(:,ii:n_obs), nhat(:,ii),entry(:,ii),s, GP );
        Kii(ii:n_obs,ii) = integral(func,0,L(ii),'ArrayValued',true)./L(ii);
        parfor_progress;
    end
else
    parfor ii = 1:n_obs
        before = ii-1;
        func = @(s) calc_dI(entry(:,ii:n_obs), exit(:,ii:n_obs), nhat(:,ii:n_obs), nhat(:,ii),entry(:,ii),s, GP );
        Kii(:,ii) = [zeros(before,1);integral(func,0,L(ii),'ArrayValued',true)./L(ii)];
        parfor_progress;
    end
end
parfor_progress(0);
Kii = Kii+tril(Kii,-1)';




end

