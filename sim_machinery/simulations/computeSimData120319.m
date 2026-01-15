function [r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData120319(R,m,uc,pnew,simtime,plotop)
if nargin<6
    plotop = 0;
end
% generate noise processes
try
    if isempty(uc) && simtime == 0
        uc = innovate_timeseries(R,m,pnew);
    end

    if simtime ~= 0
        R = setSimTime(R,simtime);
        uc = innovate_timeseries(R,m,pnew);
    end
catch
    wflag = 1;
    % disp('Cant generate noise')
end


%% Simulate New Data
% Integrate in time master fx function
try
    [xsims dum wflag] = R.IntP.intFx(R,m.x,uc,pnew,m);
catch
    disp('Simulation failed!')
    xsims{1} = nan(1,3);
    wflag = 1;
end
if wflag == 0
    try
        % Run Observer function
        % Subloop is local optimization of the observer gain
        glorg = pnew.obs.LF;
        gainlist = R.obs.glist;
        feat_sim = cell(1,length(gainlist));
        xsims_gl = cell(1,length(gainlist));
        r2mean = zeros(1,length(gainlist));
        for gl = 1:length(gainlist)
            pnew.obs.LF = glorg+gainlist(gl);
            if isfield(R.obs,'obsFx')
                [xsims_gl{gl},R,wflag(1)] = R.obs.obsFx(xsims,m,pnew,R);
            else
                xsims_gl{gl} =xsims;
            end
            if any(wflag(1))
                error('Rejection at Observation!')
            end
            % Run Data Transform d
            if isfield(R.obs,'transFx')
                [~, feat_sim{gl}, wflag(2)] = R.obs.transFx(xsims_gl{gl},R.chloc_name,R.chsim_name,1/R.IntP.dt,R.obs.SimOrd,R);
            else
                wflag(2) = 0 ;
                feat_sim{gl} = xsims_gl{gl}; % else take raw time series
            end
            % Compare Pseudodata with Real
            if isfield(R.IntP,'compFx')
                r2mean(gl)  = R.IntP.compFx(R,feat_sim{gl});
            else
                r2mean(gl) = nan;
            end
        end
        if any(wflag(2))
            error('TransFX could not compute data transform!')
        end
        [r2 ir2] = max(r2mean);
        feat_sim = feat_sim{ir2};
        xsims_gl = xsims_gl{ir2};
        pnew.obs.LF = glorg+gainlist(ir2);

        if plotop == 1
            figure(1);  R.plot.outFeatFx({R.data.feat_emp},{feat_sim},R.data.feat_xscale,R,1,[])
        end
    catch
        disp('Observation/Cost Function Failure!')
        r2 = -inf;
        ir2 =1;
        xsims_gl{1} = NaN;
        feat_sim{1} = NaN;
    end
else
    disp('Sim Output contains NaNs!')
    r2 = -inf;
    ir2 =1;
    xsims_gl{1} = NaN;
    feat_sim{1} = NaN;
end
