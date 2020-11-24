function [f] = ABC_fx_bgc_stn(x,u,P)
% state equations for a neural mass model of the basal ganglia circuit
% models the circuit between striatum, gpe, stn, gpi, and thalamus as a
% single source (no extrinsic connections)
%
% order           cells     states
% 1 = stn       - pyr       x(1,1:2)

% G(1,1) = gpe -> stn (-ve ext)

% [default] fixed parameters
%--------------------------------------------------------------------------
% G  = [2]*200;   % synaptic connection strengths
% T  = [4];               % synaptic time constants [str,gpe,stn,gpi,tha];
% R  = 2/3;                       % slope of sigmoid activation function
% NB for more pronounced state dependent transfer functions use R  = 3/2;

% input
%==========================================================================
U = u;
% time constants and intrinsic connections
%==========================================================================
T = P.T;

% intrinsic/extrinsic connections to be optimised
%--------------------------------------------------------------------------
% G = P.G; % FOR SELF CONNECTIONS 

% Motion of states: f(x)
%--------------------------------------------------------------------------
% STN: pyr

% pyramidal cells: Hidden causes - error
%--------------------------------------------------------------------------
u      =  U;
% u      =  u; %-G(:,1)*S(:,1) + u; % FOR SELF CONNECTIONS
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(1,1))./T(1,1);
% Voltage
%==========================================================================
f(:,1) = x(:,2);
f      = f'; %spm_vec(f);
