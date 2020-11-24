function R = ABCsetup_SIDoubleWell(R)
%
%% DATA SPECIFICATION
R.data.datatype{1} = 'NPD'; %%'NPD'
R.frqz = [4:.1:48];
R.frqzfull = [1:.2:48]; % legacy (only needed for SPM Structured Innovations)

R.nmsim_name = {'X1','X2','X3'}; %modules (fx) to use. These must match those listed in the fx_compile function
R.chdat_name = {'X1','X2','X3'}; % observed channels (redundant)
% R.datinds = 1:4; % Specify this when you deal with the data - ensure its not the wrong order!
R.chsim_name = {'X1','X2','X3'}; % simulated channel names
R.siminds = 1:3; % Maps from simulation to data
R.condnames = {'cond1'}; % VERIFY!!!
% Spectral characteristics
R.obs.csd.df = 0.05;
R.obs.csd.reps = 32; %96;

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @doubleWellFx;
R.IntP.compFx= @compareData_180520;

R.IntP.dt = .01;
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.obs.SimOrd = 10;
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

%% OBSERVATION
% observation function
R.obs.obsFx = @observe_data;
R.obs.gainmeth = {};
R.obs.glist =0; %
R.obs.brn =2; % 2; % burn in time

% Data Features
% fx to construct data features
R.obs.transFx = @constructGenCrossMatrix;
% These are options for transformation (NPD)
R.obs.logscale = 0;
R.obs.trans.zerobase = 0;
R.obs.trans.norm = 0;
R.obs.trans.normcat = 0;
R.obs.trans.logdetrend = 0;
R.obs.trans.gauss3 = 0;
R.obs.trans.gausSm = 0; %
R.obs.trans.interptype = 'linear';
R.obs.trans.npdscalar = 1; % optional NPD rescaling
%% OBJECTIVE FUNCTION
R.objfx.feattype = 'magnitude'; %%'ForRev'; %
R.objfx.specspec = 'auto'; %'npd'; %%'auto'; % which part of spectra to fit

%% OPTIMISATION
R.SimAn.pOptList = {'.theta'}; %,'.alist'}; %,'.flist'};
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchMax = 200;
R.SimAn.convIt.dEps = 1e-8;
R.SimAn.convIt.eqN = 5;
R.analysis.modEvi.N  = 500;
R.SimAn.scoreweight = [1 0]; %1e-8];
R.SimAn.rep = 256; % Repeats per temperature
% R.SimAn.saveout = 'xobs1';
R.SimAn.jitter = 1; % Global precision
%% PLOTTING
R.plot.outFeatFx = @genplotter_200420; 
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;


%% Perform Some Checks
if R.frqz(end)> (1/(2*R.IntP.dt))
warning('Your target frequencies are greater than Nyquist, increase your integration sim step size!')
end


