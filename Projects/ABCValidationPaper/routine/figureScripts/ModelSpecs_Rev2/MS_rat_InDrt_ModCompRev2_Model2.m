function [R p m uc] = MS_rat_InDrt_ModCompRev2_Model2(R)
%% Revised Model Space %%
% Model 1.2
%% MODEL 2 %%%
m.m = 4; % # of sources
m.x = {[0 0 0 0 0 0 0 0] [0 0]  [0 0] [0 0]}; % Initial states
m.Gint = [14 1 1 1];
m.Tint = [4 1 1 1];
m.Sint = [9 2 2 1];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[0 0 0 0 0 0 1 0] [1 0] [1 0] [1 0]};
R.obs.outstates = find([m.outstates{:}]);
for i=1:numel(R.chloc_name)
    R.obs.obsstates(i) = find(strcmp(R.chloc_name{i},R.chsim_name));
end

% Precompute xinds to make things easier with indexing
% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
        xinds(i,1) = 1;
        xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;

% setup exogenous noise
% m.uset.p = DCM.Ep;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 5e1; %.*R.InstP.dt;
uc = innovate_timeseries(R,m);


%% Prepare Priors
% 1 MMC
% 2 STR
% 3 GPE
% 4 STN

% Excitatory connections
p.A{1} =  repmat(-32,m.m,m.m);
p.A{1}(2,1) = 0; % MMC -> STR
p.A{1}(3,4) = 0; % STN -> GPE
p.A_s{1} = repmat(1/4,m.m,m.m);

p.A{2} =  repmat(-32,m.m,m.m);
p.A{2}(3,2) = 0; % STR -| GPe
p.A{2}(4,3) = 0; % GPe -| STN
p.A_s{2} = repmat(1/4,m.m,m.m);

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(1/16,size(p.C));

% Leadfield
p.obs.LF = [1 1 1 1];
p.obs.LF_s = repmat(1/4,size(p.obs.LF));

p.obs.mixing = [1]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0,size(p.obs.mixing));

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(1/16,size(p.D));

% Sigmoid transfer for connections
p.S = [0 0];
p.S_s = [1/8 1/8];

% time constants and gains
for i = 1:m.m
    if i == 1
        prec = 1/4;
    else
        prec = 1/8;
    end
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(prec,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(prec,size(p.int{i}.G));
    p.int{i}.S = zeros(1,m.Sint(i));
    p.int{i}.S_s = repmat(prec,size(p.int{i}.S));
    %     p.int{i}.BT = zeros(1,m.Tint(i));
    %     p.int{i}.BT_s = repmat(prec,size(p.int{i}.T));
end