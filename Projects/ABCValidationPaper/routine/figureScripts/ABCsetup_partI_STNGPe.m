function R = ABCsetup_partI_STNGPe(R)
%
%% DATA SPECIFICATION
R.data.datatype{1} = 'NPD'; % data type
R.frqz = [6:.2:48]; % frequency support
R.frqzfull = [1:.2:200]; % legacy (only needed if using SPM Structured Innovations)

R.nmsim_name = {'GPE','STN'}; %modules (fx) to use. These must match those listed in the fx_compile function
R.chdat_name = {'GPE','STN'}; % observed channels (redundant)
R.chsim_name = {'GPE','STN'}; % simulated channel names
R.siminds = 1:2; % Maps from simulation to data
R.condnames = {'OFF'}; % 

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @ABC_fx_compile_120319; % main dynamics function

R.IntP.dt = .001; % integration step size
R.IntP.Utype = 'white_covar'; % innovation type
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays


%% OBSERVATION
R.obs.obsFx = @observe_data; % observation function
R.obs.gainmeth = {'obsnoise','unitvar'}; % options for observation model
R.obs.glist =0; % gain sweep optimization range [min max listn] (log scaling)
R.obs.brn =2; % burn in time
LF = [1 1]*10; % leadfield
R.obs.LF = LF;
R.obs.Cnoise = [0.2 0.2]; % Sensor Noise SNR prior i.e. 1/x signal to noise ratio

% Spectral characteristics (determine sim length)
R.obs.csd.df = 0.5;
R.obs.csd.reps = 32;
R = setSimTime(R,R.obs.csd.reps);

%% Construct Data Features
R.obs.transFx = @constructGenCrossMatrix; % fx to construct data features

% These are options for transformation (NPD)
R.obs.logscale = 0;
R.obs.trans.zerobase = 0;
R.obs.trans.norm = 0;
R.obs.trans.normcat = 0;
R.obs.trans.logdetrend = 0;
R.obs.trans.gauss3 = 0;
R.obs.trans.gausSm = 4; %  Gaussian smoothing 
R.obs.trans.interptype = 'linear';
R.obs.trans.npdscalar = 1; % optional NPD rescaling

%% OBJECTIVE FUNCTION
R.objfx.feattype = 'magnitude'; %%'ForRev'; %
R.objfx.specspec = 'npd'; %'npd'; %%'auto'; % which part of spectra to fit
R.objfx.compFx= @compareData_180520; % cost function

%% OPTIMISATION
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A','.D','.obs.Cnoise'}; % Free parameters
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2); % only used for constructing support for plotting
R.SimAn.searchMax = 200; % maximum number of iterations
R.SimAn.convIt.dEps = 5e-2; % change of objective fx at convergence
R.SimAn.convIt.eqN = 5; % maximum number of attempts
R.SimAn.scoreweight = [1 0]; % complexity weighted goodness of fit [acc cmplx]
R.SimAn.rep = 512; % Draws per iteration
R.SimAn.jitter = 1; % Global precision

%% ANALYSIS
R.analysis.modEvi.N  = 500; % Number of draws

%% PLOTTING
R.plot.outFeatFx = @genplotter_200420; 
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;






