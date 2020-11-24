function [f] = ABC_fx_bgc_gpi(x,u,P)
% state equations for a neural mass model of the basal ganglia circuit
% models the circuit between striatum, gpe, stn, gpi, and thalamus as a
% single source (no extrinsic connections)
%
% order           cells     states
% 1 = gpi       - ii        x(1,7:8)

% G(1,1) = str -> gpi (-ve ext)
% G(1,2) = stn -> gpi (+ve ext)
% G(1,3) = gpi -> gpi (-ve self)

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = P.Rz(2:end);              % gain of activation function (1st is extrinsic- so remove)
S = sigmoidin(x,R,0);
S = S';

% R    = R.*exp(P.S);              % gain of activation function
% F    = 1./(1 + exp(-R*x + 0));   % firing rate
% S    = F - 1/(1 + exp(0));       % deviation from baseline firing (0)

% input
%==========================================================================
U = u;

% time constants and intrinsic connections
%==========================================================================
T = P.T;

% intrinsic/extrinsic connections to be optimised
%--------------------------------------------------------------------------
G = P.G;

% Motion of states: f(x)
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 - GPi: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  U;
% u      =  -G(:,1)*S(:,1) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(1,1))./T(1,1);

% Voltage
%==========================================================================
f(:,1) = x(:,2);
f      = f'; %spm_vec(f);




