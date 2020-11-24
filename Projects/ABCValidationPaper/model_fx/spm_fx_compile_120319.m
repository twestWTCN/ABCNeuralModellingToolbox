function [xstore_cond,tvec,wflag,J,Es] = spm_fx_compile_120319(R,x,uc,pc,m)
% To Do:
% 1)Precompute the expectations of the within source parameters and take
%   outside of the integration loop.
% If you want to estimate the noise floor - then use 'decon' to deconnect
% both intrinsic/extrinsic couplings.
if isfield(R.IntP,'getNoise') && R.IntP.getNoise == 1
    decon = 0;
else
    decon = 1;
end

cs = 0; % cond counter
wflag= 0; tvec = [];
for condsel = 1:numel(R.condnames)
    cs = cs+1;
    us = uc{cs};
    p = pc;
    % Compiles NMM functions with delays, seperates intrinsic and extrinsic
    % dynamics then summates
    xinds = m.xinds;
    % model-specific parameters
    %==========================================================================
    % model or node-specific state equations of motions
    %--------------------------------------------------------------------------
    fx{1} = @spm_fx_erp;                                    % ERP model
    fx{2} = @spm_fx_cmc_local;                                    % CMC model
    fx{3} = @spm_fx_bgc;                                    % basal ganglia circuit
    %     fx{4} = @spm_fx_mmc;                                    % motor micro circuit (SPM version)
    
    fx{4} = @ABC_fx_bgc_mmc;                                    % motor micro circuit
    fx{5} = @ABC_fx_bgc_str;
    fx{6} = @ABC_fx_bgc_gpe;
    fx{7} = @ABC_fx_bgc_stn;
    fx{8} = @ABC_fx_bgc_gpi;
    fx{9} = @ABC_fx_bgc_thal;
    % indices of extrinsically coupled hidden states
    %--------------------------------------------------------------------------
    efferent(1,:) = [9 9 9 9];               % sources of ERP connections
    afferent(1,:) = [4 8 5 8];               % targets of ERP connections
    
    efferent(2,:) = [3 3 7 7];               % sources of CMC connections
    afferent(2,:) = [2 8 4 6];               % targets of CMC connections
    
    efferent(3,:) = [9 9 9 9];               % sources of BGC connections (thalamus)
    afferent(3,:) = [2 6 2 6];               % targets of BGC connections (striatum & STN)
    
    efferent(4,:) = [3 3 6 7];               % ORIG sources of MMC connections
    afferent(4,:) = [2 4 8 0];               % targets of MMC connections
    %     efferent(4,:) = [7 7 7 7];               % sources of MMC connections
    %     afferent(4,:) = [8 8 8 8];                  % forward deep/middle; back deep/superficial
    efferent(5,:) = [1 1 1 1];               % sources of STR connections
    afferent(5,:) = [2 2 2 2];               % targets of STR connections
    
    efferent(6,:) = [1 1 1 1];               % sources of GPE connections
    afferent(6,:) = [2 2 2 2];               % targets of GPE connections
    
    efferent(7,:) = [1 1 1 1];               % sources of STN connections
    afferent(7,:) = [2 2 2 2];               % targets of STN connections
    
    efferent(8,:) = [1 1 1 1];               % sources of GPI connections
    afferent(8,:) = [2 2 2 2];               % targets of GPI connections
    
    efferent(9,:) = [1 1 1 1];               % sources of THAL connections
    afferent(9,:) = [2 2 2 2];               % targets of THAL connections
       
    
    % scaling of afferent extrinsic connectivity (Hz)
    %--------------------------------------------------------------------------
    E(1,:) = [1 0 1 0]*200;                    % ERP connections
    E(2,:) = [1 .3571 1 .625]*100000;          % CMC connections (to ctx) with T = [2 2 16 28] gives [200 100 200 100] = regular DCM
    E(3,:) = [1.8 1.2 1.8 1.2]*10000;         % BGC connections (to bgc) with T_str=8 and T_stn=4 gives A = 144 and 48
    %     E(4,:) = [.2 .2 -.2 -.2]*(200);  %500           % MMC connections (to mmc) with T_mp=3 and T_sp=2 gives A = 270 and 180; with T_dp=18 gives A=200
    E(4,:) = [.2 .2 -.2 -.2]*10000;
    %% to calculate E divide the target value for A by the value of the time constant (in seconds, i.e. 0.018)
    % E(5,:) = [.5 .5 -.5 -.5]*100000;             % STR connections
    % E(6,:) = [.5 .5 -.5 -.5]*100000;             % GPE connections
    % E(7,:) = [ 1  1 -.1  -1]*100000;             % STN connections
    % E(8,:) = [.5 .5 -.5 -.5]*100000;               % GPI connections
    % E(9,:) = [.5 .5 -.5 -.5]*100000;               % THAL connections
    
    E(5,:) = [.2 .2 -.2 -.2]*8000;             % STR connections
    E(6,:) = [.2 .2 -.2 -.2]*10000;             % GPE connections
    E(7,:) = [.2 .2 -.2 -.2]*10000;             % STN connections
    E(8,:) = [.2 .2 -.2 -.2]*8000;             % GPI connections
    E(9,:) = [.2 .2 -.2 -.2]*5000;  %500       % THAL connections
    
    % get the neural mass models {'ERP','CMC'}
    %--------------------------------------------------------------------------
    n     = m.m;
    model = m.dipfit.model;
    for i = 1:n
        if  strcmp(model(i).source,'ERP')
            nmm(i) = 1;
        elseif strcmp(model(i).source,'CMC')
            nmm(i) = 2;
        elseif strcmp(model(i).source,'BGC')
            nmm(i) = 3;
        elseif strcmp(model(i).source,'MMC')
            nmm(i) = 4;
        elseif  strcmp(model(i).source,'STR')
            nmm(i) = 5;
        elseif  strcmp(model(i).source,'GPE')
            nmm(i) = 6;
        elseif  strcmp(model(i).source,'STN')
            nmm(i) = 7;
        elseif  strcmp(model(i).source,'GPI')
            nmm(i) = 8;
        elseif  strcmp(model(i).source,'THAL')
            nmm(i) = 9;
        end
    end
    
    %% Pre-integration extrinsic connection parameters
    
    % Compute value of delays from lognormal mean
    D = zeros(m.m);
    D(p.D>-30) = 4/1000; % set all delay priors to 4ms.
    
    D(2,1) = 3/1000;   % M1 to STR (Jaeger and Kita, 2011)
    D(4,1) = 3/1000;  % M1 to STN (Jaeger and Kita, 2011)
   
    D(3,2) = 7/1000;   % STR to GPe (Kita and Kitai 1991)
    D(5,2) = 12/1000;  % STR to GPi (Kita and Kitai 1991)

    D(4,3) = 1/1000;    % GPe to STN (Jaeger and Kita, 2011)
    D(5,3) = 1/1000;    % GPe to GPi (Jaeger and Kita, 2011)
   
    D(3,4) = 3/1000;    % STN to GPe (Kita and Kitai 1991)
    D(5,4) = 3/1000;    % STN to GPi (Kita and Kitai 1991)
    
    D(6,5) = 3/1000;    % GPi to Thal (Stoelzel J Neurosci. 2017)
    
    D(1,6) = 3/1000;   % Thal to M1 (Lumer, Edelman, Tononi; 1997)
    D(1,6) = 8/1000;   % M1 to Thal (Lumer, Edelman, Tononi; 1997)
    
    
    D = D(1:m.m,1:m.m);
    D = ceil(D.*exp(p.D).*(1/R.IntP.dt)); % As expectation of priors and convert units to steps
    D(D<((1e-3)/R.IntP.dt)&D>0) = floor((2e-3)/R.IntP.dt); % Minimum 1ms
    
    if (R.IntP.buffer-max(max(D)))<=0
        R.IntP.buffer = max(max(D)) + 2;
        disp(['Delay is bigger than buffer, increasing buffer to: ' num2str(R.IntP.buffer)])
    end
    if R.IntP.buffer > 1e3
        disp('Delays are implausibly large!')
        wflag = 1;
        break
    end
    
    
    Ds = zeros(size(D));Dt = zeros(size(D));
    % Now find indices of inputs
    % Currently no seperation between inh and excitatory
    for i = 1:length(nmm) % target
        for j = 1:length(D(i,:)) % source
            if D(i,j)>0
                Ds(i,j) = efferent(nmm(j),1); % input sources
                Ds(i,j) = (m.xinds(j,1)-1)+Ds(i,j);
                Dt(i,j) = afferent(nmm(i),1); % input targets
                Dt(i,j) = (m.xinds(j,1)-1)+Dt(i,j);
            end
        end
    end
       
    % Condition Dependent Modulation of Synaptic gains
    %-----------------------------------------
    for i = 1:m.m
        if cs ~= R.Bcond
            p.int{i}.T = p.int{i}.T;
        else
            p.int{i}.T = p.int{i}.T + p.int{i}.BT;
        end
    end
    
    % Rescale background Input
    for i = 1:m.m
        C    = exp(p.C(i));
        us(:,i) = C*us(:,i)*0.01;
    end
    % Extrinsic connections
    %--------------------------------------------------------------------------
    %
    % alist = [1 2; 3 4];
    alist = [1; 3];
    for i = 1:numel(p.A)
        if cs ~= R.Bcond
            A{i} = decon*exp(p.A{i});
        else
            A{i} = decon*exp(p.A{i}+p.B{i}); % Add the second condition
        end
        %     A{alist(i,2)} = exp(p.A{i});
    end
    
    % detect and reduce the strength of reciprocal (lateral) connections
    %--------------------------------------------------------------------------
    % TOL   = exp(2);
    % for i = 1:numel(A)
    %     L    = (A{i} > TOL) & (A{i}' > TOL);
    %     A{i} = A{i}./(1 + 4*L);
    % end
    
    % and scale of extrinsic connectivity (Hz)
    %--------------------------------------------------------------------------
    for j = 1:n
        for i = 1:n
            for k = 1:numel(p.A)
                A{k}(i,j) = E(nmm(i),alist(k,1))*A{k}(i,j);
            end
        end
    end
    
    % synaptic activation function priors
    %--------------------------------------------------------------------------
    Rz_base     = 2/3;                      % gain of sigmoid activation function   
    B = 0;
    %% Precompute parameter expectations
    % Parameter Priors
    pQ = getModelPriors(m);
    
    nbank = cell(1,n); qbank = cell(1,n);
    for i = 1:n
        N.x  = m.x{i};
        nbank{i} = N;
        Q = p.int{i};
        Q.T = pQ(i).T.*exp(Q.T);
        Q.G = decon*pQ(i).G.*exp(Q.G);
        Q.Rz = Rz_base.*exp(Q.S);
        Rz(i) = Q.Rz(1);
        Q.C  = p.C(i,:);
        qbank{i} = Q;
    end
    
    %% TIME INTEGRATION STARTS HERE ===========================================
    f = zeros(xinds(end),1); dt = R.IntP.dt;
    if iscell(x)
        xstore= full(repmat(spm_vec(x),1,R.IntP.buffer));
    else
        xstore = x;
    end
    % pad out the rest of xstore with zeros
    %     xstore =    [xstore zeros(m.xinds(end),(R.IntP.nt+1)-size(xstore,2))];
    
    %     xstore = [xstore nan(size(xstore,1),R.IntP.nt-R.IntP.buffer)];
    xint = zeros(m.n,1);
    TOL = exp(-4);
    for tstep = R.IntP.buffer:R.IntP.nt
        % assemble flow
        %==========================================================================
        N     = m;
        for i = 1:n % targets
            fA = [];
            % extrinsic flow
            %----------------------------------------------------------------------
            for j = 1:n % sources
                for k = 1:numel(p.A) % connection type
                    if abs(A{k}(i,j)) > TOL
                        %                 ik       = afferent(nmm(i),k);
                        %                 jk       = efferent(nmm(j),k);
                        %                 xd = spm_unvec(xback(:,end-D(i,j)),M.x);
                        xD = xstore(Ds(i,j),tstep-D(i,j));
                        %                     fA(Dt(i,j)) = A{k}(i,j)*S(xD,Rz,B);
                        fA = [fA  A{k}(i,j)*sigmoidin(xD,Rz(j),B)]; % 1st Rz is slope!
                    end
                end
            end
            % intrinsic flow at target
            %----------------------------------------------------------------------
            ui   = us(tstep,i)+ sum(fA);
            xi = xstore(m.xinds(i,1):m.xinds(i,2),tstep)';
            f(m.xinds(i,1):m.xinds(i,2)) = fx{nmm(i)}(xi,ui,qbank{i});
            %             f(Dt(1,i))  = f(Dt(1,i)) + sum(fA) ;
        end
        xint = xint + (f.*dt);
        xstore = [xstore xint]; % This is done for speed reasons! Faster than indexing (!!)
        if any(isnan(xint))
            a = 1;
        end
        if tstep >R.IntP.buffer*10
            if any(xint>1e4) || any(isnan(xint))
                wflag= 1;
                break
            end
            pp1 = 1;
        end
        % disp(tstep/R.IntP.nt)
        % xint= spm_unvec(x,M.x);
    end
    if wflag == 1
        xstore_cond{condsel} = NaN;
    end
    xstore_cond{condsel} = xstore;
    if nargout>3
        [J{condsel},Es{condsel}] = findJacobian(R,xstore(:,end-R.IntP.buffer:end),uc,p,m);
    end    % tvec = linspace(R.IntP.buffer*R.IntP.dt,R.IntP.nt*R.IntP.dt,R.IntP.nt);
    a = 1;
end

