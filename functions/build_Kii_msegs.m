function [ Kii ] = build_Kii_msegs( entry, exit, n_obs, GP,nsegs )

% code is only set up to work for 2 segs maximum atm
% segment 1 with segment 1

d = exit(1:2,:)-entry(1:2,:);
L1 = hypot(d(1,:),d(2,:));
[ Kii ] = build_Kii( entry(1:2,:), exit(1:2,:), n_obs, GP ).*bsxfun(@times,L1',L1);



% segment 2 with segment 2
Isegs = (nsegs >= 2);  
n_ss = sum(Isegs);
ss = 2;
d = exit(3:4,Isegs)-entry(3:4,Isegs);
L2 = hypot(d(1,:),d(2,:));
Kii(Isegs,Isegs) = Kii(Isegs,Isegs) + build_Kii( entry(2*ss-1:2*ss,Isegs),...
    exit(2*ss-1:2*ss,Isegs), n_ss, GP ).*bsxfun(@times,L2',L2);

% segment 2 with segment 1 
% ss = 2
entry12 = [entry(1:2,:) entry(3:4,Isegs)];
exit12 = [exit(1:2,:) exit(3:4,Isegs)];

d = exit12-entry12;
L12 = hypot(d(1,:),d(2,:));
tmp = build_Kii(entry12,exit12,length(d),GP).*bsxfun(@times,L12',L12);
Kii(:,Isegs) = Kii(:,Isegs) + tmp(1:length(L1),length(L1)+1:end);
Kii(Isegs,:) = Kii(Isegs,:) + tmp(length(L1)+1:end,1:length(L1));


% correct for total length
LT = L1; LT(Isegs) = LT(Isegs) + L2;
Kii = Kii./bsxfun(@times,LT',LT);

%% TODO correct for the lengths

end

