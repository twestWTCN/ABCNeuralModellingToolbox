function R = ABCsetup_partII_FullModel(R)
%
%% DATA SPECIFICATION
R.data.datatype{1} = 'NPD'; %%'NPD'
R.frqz = [6:.2:48];
R.frqzfull = [1:.2:200]; % legacy (only needed for SPM Structured Innovations)

R.nmsim_name = {'MMC','STR','GPE','STN'}; %modules (fx) to use. These must match those listed in the fx_compile function
R.chdat_name = {'MMC','STR','GPE','STN'}; % observed channels (redundant)
% R.datinds = 1:4; % Specify this when you deal with the data - ensure its not the wrong order!
R.chsim_name = {'MMC','STR','GPE','STN'}; % simulated channel names
R.siminds = 1:4; % Maps from simulation to data
R.condnames = {'OFF'}; % VERIFY!!!
% Spectral characteristics
R.obs.csd.df = 0.5;
R.obs.csd.reps = 32; %96;

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @spm_fx_compile_120319;
R.IntP.compFx= @compareData_180520;

R.IntP.dt = .001;
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.obs.SimOrd = 8;
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

%% OBSERVATION
% observation function
R.obs.obsFx = @observe_data;
R.obs.gainmeth = {'obsnoise','unitvar'}; %,'submixing'}; %,'lowpass'}; ,'leadfield' %unitvar'mixing'
R.obs.glist =0; %linspace(-5,5,12);  % gain sweep optimization range [min max listn] (log scaling)
R.obs.brn = 2; % 2; % burn in time
LF = [1 1 1 1]*10; % Fit visually and for normalised data
R.obs.LF = LF;
R.obs.Cnoise = [0.2 0.2 0.2 0.2]; % Noise gain on the observation function

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
R.obs.trans.gausSm = 4; % This is off but is switched on to 1 Hz at data processing stage
R.obs.trans.interptype = 'linear';
R.obs.trans.npdscalar = 1; % optional NPD rescaling
%% OBJECTIVE FUNCTION
R.objfx.feattype = 'magnitude'; %%'ForRev'; %
R.objfx.specspec = 'npd'; %%'auto'; % which part of spectra to fit

%% OPTIMISATION
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A','.D'}; %,'.S','.int{src}.G','.int{src}.S','.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchMax = 200;
R.SimAn.convIt.dEps = 1e-6;
R.SimAn.convIt.eqN = 5;
R.analysis.modEvi.N  = 500;
R.SimAn.scoreweight = [1 1e-6];
R.SimAn.rep = 256; %512; % Repeats per temperature
% R.SimAn.saveout = 'xobs1';
R.SimAn.jitter = 1; % Global precision
%% PLOTTING
R.plot.outFeatFx = @genplotter_200420; 
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;






