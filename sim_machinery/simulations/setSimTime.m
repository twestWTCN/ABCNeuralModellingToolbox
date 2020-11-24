function R = setSimTime(R,N)
disp('You are reinitializing the simulation time!')
disp(sprintf('The Simulation was %.2f seconds long',R.IntP.tend))

% R.obs.csd.df = 0.5;
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = (log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
% R.IntP.dt = .0005;
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = floor(R.IntP.tend/R.IntP.dt);
dfact = fsamp/(2*2^(R.obs.SimOrd));
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

tvec_obs = R.IntP.tvec;
tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
R.IntP.tvec_obs = tvec_obs;

disp(sprintf('The Simulation is now %.2f seconds long',R.IntP.tend))

disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));
