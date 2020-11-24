function R = prepareDoubleWellData(R,plotop,nulldat)
if nargin<2
    plotop = 1;
end
if nargin<3
    nulldat = 0;
end

if nulldat == 1
    disp('Prepping Null Data')
end
% prepareratdata_group(R.rootn,R.projectn);
load([R.path.datapath  '\DoubleWell\sim1_f14_30_theta_3_neg2_1dot5.mat']);
% x = (x-mean(x))./std(x);
[F feat_out] = R.obs.transFx({x'},1:3,1/dt,R.obs.SimOrd ,R); 

% Set data as working version
R.data.feat_emp = feat_out;
% squeeze(meannpd_data(1,1,1,1,:))
R.data.feat_xscale{1} = F{1};

% Plot CSD
if plotop ==1
     R.plot.outFeatFx({R.data.feat_emp},{R.data.feat_emp},R.data.feat_xscale,R,1,[])
end


