function [Kii] = Kii_integral1_segment( entry, exit,nhat, n_obs,nSegs, GP,parallel )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

L = hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));
% L2 = hypot(exit(3,:)-entry(3,:),exit(4,:)-entry(4,:));
% Lt = L;
% Lt() = 
% nhat = (exit-entry)./L;
Kii = zeros(n_obs,n_obs);

parfor_progress(n_obs);
if ~parallel  
    for ii = 1:n_obs
        Lt = 0;
        for tt = 1:nSegs(ii)
            L = hypot(exit(2*tt-1,ii)-entry(2*tt-1,ii),exit(2*tt,ii)-entry(2*tt,ii));
            Lt = Lt +L;
            func = @(s) calc_dI_segs(entry(:,ii:n_obs), exit(:,ii:n_obs), nhat(:,ii:n_obs),...
                nhat(:,ii),entry(2*tt-1:2*tt,ii),s,nSegs(ii:n_obs), GP );
            Kii(ii:n_obs,ii) =Kii(ii:n_obs,ii)+ integral(func,0,L,'ArrayValued',true);
        end
        Kii(ii:n_obs,ii) = Kii(ii:n_obs,ii)./Lt;
        parfor_progress;
    end
else
    parfor ii = 1:n_obs
        Lt = 0;
        for tt = 1:nSegs(ii)
            before = ii-1;
            L = hypot(exit(2*tt-1,ii)-entry(2*tt-1,ii),exit(2*tt,ii)-entry(2*tt,ii));
            Lt = Lt +L;
            func = @(s) calc_dI_segs(entry(:,ii:n_obs), exit(:,ii:n_obs), nhat(:,ii:n_obs),...
                nhat(:,ii),entry(2*tt-1:2*tt,ii),s,nSegs(ii:n_obs), GP );
            Kii(:,ii) =Kii(:,ii)+ [zeros(before,1);integral(func,0,L,'ArrayValued',true)];
        end
        Kii(:,ii) = Kii(:,ii)./Lt;
        parfor_progress;
    end
end
parfor_progress(0);
Kii = Kii+tril(Kii,-1)';




end

