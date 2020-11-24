function    S = sigmoidin(x,Rz,B)
S     = (1./(1 + exp(-Rz(:).*x(:) + B))) - 1/(1 + exp(B));