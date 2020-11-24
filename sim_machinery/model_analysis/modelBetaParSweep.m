function parsweep = modelBetaParSweep(m,p,parsweep,R)
% R.obs.gainmeth = R.obs.gainmeth(1:2);
% R.obs.gainmeth = {'unitvar','leadfield','submixing'};
bwid = 1;
% R.obs.csd.df = 1;
% R.obs.csd.reps = 96; %32;

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;

u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);


Qlist = parsweep.Qlist;
Rlist = parsweep.Rlist;
% Put inside loop to reduce memory load
% for q = 1:length(parsweep.Qlist)
%     for r = 1:length(parsweep.Rlist)
%         pa = p;
%         eval(['pa' parsweep.Q ' = Qlist(q);'])
%         eval(['pa' parsweep.R ' = Rlist(r);'])
%         psweep{q,r} = pa;
%     end
% end
betaPowBank = nan(m.m,length(parsweep.Qlist),length(parsweep.Rlist));
maxfrqBank = nan(m.m,length(parsweep.Qlist),length(parsweep.Rlist));
maxPowBank = nan(m.m,length(parsweep.Qlist),length(parsweep.Rlist));
% frqPowBank = cell(length(parsweep.Qlist),length(parsweep.Rlist));
for q = 1:length(parsweep.Qlist)
    parfor r = 1:length(parsweep.Rlist)
        frqlist = 8:0.5:48;
        pnew =p;
        pnew = pareval(['subject' parsweep.Q ' = input;'],pnew,parsweep.Q,Qlist(q));
        pnew = pareval(['subject' parsweep.R ' = input;'],pnew,parsweep.R,Rlist(r));
        % Integrate
        xsims = R.IntP.intFx(R,m.x,u,pnew,m);
        % Run Observer function
        if isfield(R.obs,'obsFx')
            [xsims,~,wflag] = R.obs.obsFx(xsims,m,pnew,R);
        end
        if wflag == 0
            % Plot
            %         figure
            %         plotRawTimeSeries(R,xsims{1})
            %         set(gcf,'Position',[680  281  824  696])
            frqPow = zeros(m.m,length(frqlist));
            for i = 1:m.m
                for j = 1:length(frqlist)
                    frqPow(i,j) = bandpower(xsims{1}(i,:),1/R.IntP.dt,[frqlist(j)-bwid frqlist(j)+bwid]);
                end
            end
            betaPow = sum(frqPow(:,frqlist>=14 & frqlist<22),2); % beta-band (lower) power
            [maxpow loc] = max(frqPow,[],2);
            maxfrqBank(:,q,r) = [frqlist(loc(1)) frqlist(loc(2))]; % frequency of max power
            maxPowBank(:,q,r) = [maxpow(1) maxpow(2)]; % Entire frq bank
            betaPowBank(:,q,r) = betaPow;
        end
        disp([q r])
    end
end

parsweep.Rlist = Rlist;
parsweep.Qlist = Qlist;
% parsweep.R = '.A{1}(3,4)';
% parsweep.Q = '.A{1}(4,3)';
parsweep.frqPowBank = maxPowBank;
parsweep.betaPowBank = betaPowBank;
parsweep.maxfrqBank = maxfrqBank;

