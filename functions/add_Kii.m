function [ Kii_new ] = add_Kii( entry_new, exit_new, entry_old, exit_old, Kii_old, GP,options )
%[ Kii_new ] = add_Kii( entry_new, exit_new, entry_old, exit_old, Kii_old, GP )
%   Detailed explanation goes here

if nargin < 7
    options = 'on';
end
if strcmp(options,'on')
    disp('Adding observations to Kii')
end

[~,n_old] = size(entry_old);
[~,n_new] = size(entry_new);
n_obs = n_old + n_new;

sig_f1 = GP.sig_f1;
M = GP.M;
a = GP.a;
b = GP.b;

obs = [entry_old(1,:), entry_new(1,:);
    exit_old(1,:), exit_new(1,:);
    entry_old(2,:), entry_new(2,:);
    exit_old(2,:), exit_new(2,:)];

%% Line-integral/line-integral covariance
% parallel computation
obs1 = zeros(3,n_obs*((n_obs-1)/2+1));
obs2 = obs1;

for ww=1:4
    aa = repmat(obs(ww,:),n_obs,1);
    obs1(ww,:) = aa(logical(triu(ones(n_obs))))';
%     obs1(ww,:) = aa(:,end)';
    aa= aa'; 
    obs2(ww,:) = aa(logical(triu(ones(n_obs))))';
%     obs2(ww,:) = aa(:,end)';
end
obs11 = obs1([1 3],n_old*((n_old-1)/2+1)+1:end); 
obs12 = obs1([2 4],n_old*((n_old-1)/2+1)+1:end);
obs21 = obs2([1 3],n_old*((n_old-1)/2+1)+1:end); 
obs22 = obs2([2 4],n_old*((n_old-1)/2+1)+1:end);

%% Gram matrix components
n_add = n_obs*((n_obs-1)/2+1)-n_old*((n_old-1)/2+1);
K_ii_new = zeros(n_add,1);
if strcmp(options,'on')
    parfor_progress(n_obs*((n_obs-1)/2+1));
end
parfor qq=1:n_add % utilize loop independence

    Ltot1 = norm(obs12(:,qq)-obs11(:,qq));  % length of ray 1
    Ltot2 = norm(obs22(:,qq)-obs21(:,qq));  % length of ray 2
    % obs11 is the entry point of ray 1, obs12 is the exit point of ray 1
    % obs21 is the entry point of ray 2, obs22 is the exit point of ray 2
    
    n1=(obs12(:,qq)-obs11(:,qq))/norm(obs12(:,qq)-obs11(:,qq)); % nhat of ray 1
    n2=(obs22(:,qq)-obs21(:,qq))/norm(obs22(:,qq)-obs21(:,qq)); % nhat of ray 2
    nx1=n1(1); ny1=n1(2); nx2=n2(1); ny2=n2(2);     
    n1_hat_1 = M.*n1;
    n2_hat_1 = M.*n2;
    v = obs11(:,qq)-obs21(:,qq);
    v_hat_1  = M.*v;
    x10 = obs11(:,qq); x11 = obs12(:,qq);   % entry and exit points of line 1
    x20 = obs21(:,qq); x21 = obs22(:,qq);   % entry and exit points of line 2

    L1s = 0;                    % start location of the integral
    L1  = norm(x11-x10);        % finish location of the integral

    L2s = 0;                    % start location of second integral
    L2  = norm(x21-x20);        % finish location of second integral


    %% brutal alternative
    k1 = @(s,t) sig_f1^2*exp(-0.5*( n1'*n1_hat_1*s.^2 - 2*n1'*n2_hat_1*t.*s +n2'*n2_hat_1*t.^2 +2*v'*n1_hat_1*s-2*v'*n2_hat_1.*t + v'*v_hat_1 ));

    k11 = @(s,t) ( (( s*nx1+x10(1)-(t*nx2+x20(1)) ).^4*M(1)^2-6*M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2+3)*M(1)^2 -...
        (a/b)*2*(1-M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2).*(1-M(2)*(s*ny1+x10(2)-(t*ny2+x20(2))).^2)*M(1)*M(2) +...
        (a/b)^2*(( s*ny1+x10(2)-(t*ny2+x20(2)) ).^4*M(2)^2-6*M(2)*( s*ny1+x10(2)-(t*ny2+x20(2)) ).^2+3)*M(2)^2 ).*k1(s,t);

    k22 = @(s,t) ((a+b)^2/(b^2))*(1-M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2).*(1-M(2)*(s*ny1+x10(2)-(t*ny2+x20(2))).^2)*M(1)*M(2).*k1(s,t);

    k33 = @(s,t) ( (a/b)^2*(( s*nx1+x10(1)-(t*nx2+x20(1)) ).^4*M(1)^2-6*M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2+3)*M(1)^2 -...
        (a/b)*2*(1-M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2).*(1-M(2)*(s*ny1+x10(2)-(t*ny2+x20(2))).^2)*M(1)*M(2) +...
        (( s*ny1+x10(2)-(t*ny2+x20(2)) ).^4*M(2)^2-6*M(2)*( s*ny1+x10(2)-(t*ny2+x20(2)) ).^2+3)*M(2)^2 ).*k1(s,t);

    k12 = @(s,t) ( ((a+b)/b)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).*(s*ny1+x10(2)-(t*ny2+x20(2))).*(( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2*M(1)-3)*M(2)*M(1)^2 -...
        (a*(a+b)/b^2)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).*(s*ny1+x10(2)-(t*ny2+x20(2))).*(( s*ny1+x10(2)-(t*ny2+x20(2)) ).^2*M(2)-3)*M(1)*M(2)^2 ).*k1(s,t);

    k23 = @(s,t) ( ((a+b)/b)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).*(s*ny1+x10(2)-(t*ny2+x20(2))).*(( s*ny1+x10(2)-(t*ny2+x20(2)) ).^2*M(2)-3)*M(1)*M(2)^2 -...
        (a*(a+b)/b^2)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).*(s*ny1+x10(2)-(t*ny2+x20(2))).*(( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2*M(1)-3)*M(2)*M(1)^2 ).*k1(s,t);

    k13 = @(s,t) ( -(a/b)*(( s*nx1+x10(1)-(t*nx2+x20(1)) ).^4*M(1)^2-6*M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2+3)*M(1)^2 +...
        ((a^2+b^2)/(b^2))*(1-M(1)*( s*nx1+x10(1)-(t*nx2+x20(1)) ).^2).*(1-M(2)*(s*ny1+x10(2)-(t*ny2+x20(2))).^2)*M(1)*M(2)-...
        (a/b)*(( s*ny1+x10(2)-(t*ny2+x20(2)) ).^4*M(2)^2-6*M(2)*( s*ny1+x10(2)-(t*ny2+x20(2)) ).^2+3)*M(2)^2).*k1(s,t);

    K_ii_new(qq) = K_ii_new(qq) + (1/(Ltot1*Ltot2))*integral2( @(s,t) ...
        nx1^2*nx2^2*k11(s,t)+4*nx1*nx2*ny1*ny2*k22(s,t)+...
        ny1^2*ny2^2*k33(s,t)+2*( nx2^2*nx1*ny1+nx1^2*nx2*ny2 )*k12(s,t)+...
        2*( ny2^2*nx1*ny1+ny1^2*nx2*ny2 )*k23(s,t)+( ny1^2*nx2^2 + nx1^2*ny2^2 )*k13(s,t)...
        ,L1s,L1s+L1,L2s,L2s+L2);
    if strcmp(options,'on')
        parfor_progress;
    end
end
if strcmp(options,'on')
    parfor_progress(0);
end
%%
Kq = Kii_old(logical(triu(ones(n_old))));
Kq_t = [Kq; K_ii_new];

Kn      = zeros(n_obs);
Kn(logical(triu(ones(n_obs)))) = Kq_t;  % put K_li into the upper triangle
Kii_new    = Kn+triu(Kn,+1)';              % add in the lower triangle


end

