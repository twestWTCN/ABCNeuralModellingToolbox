function [R p m uc] = MS_doubleWellSpec_1(R)
m.m = 3;
R.obs.outstates = 1:3
% setup exogenous noise
% m.uset.p = DCM.Ep;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 1; %.*R.InstP.dt;
uc = innovate_timeseries(R,m);

m.x = [0 0 0];
%% Prepare Priors
p.alist = [0 0];
p.alist_s = [1/4 1/4];

p.flist = [0 0];
p.flist_s = [1/8 1/8];

p.theta = [0 0 0];
p.theta_s = [1/4 1/4 1/4];

p.obs.LF = 0;