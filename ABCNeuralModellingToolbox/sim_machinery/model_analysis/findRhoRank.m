function [RhoRank RhoRed] = findRhoRank(Rho,eps);
erh = eig(Rho);
erh = erh(end:-1:1);

erh = erh./sum(erh);
erhCS = cumsum(erh);

RhoRank  = sum(erhCS<eps);
RhoRed = ((size(Rho,1)-RhoRank)./size(Rho,1))*100;