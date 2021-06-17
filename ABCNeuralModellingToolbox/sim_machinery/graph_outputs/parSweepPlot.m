function parSweepPlot(R,parsweep,cmap)
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery\graph_outputs\cmap_ABC.mat')
% cmap = linspecer(256);
cmap = brewermap(256,'RdYlBu');
cmap = cmap(end:-1:1,:);
colormap(cmap)
set(gcf,'color','w');

Xr = min(parsweep.Rlist):.01:max(parsweep.Rlist);
Xq = min(parsweep.Qlist):.01:max(parsweep.Qlist);
[X,Y] = meshgrid(Xr,Xq);
% V = interp2(parsweep.Rlist,parsweep.Qlist,squeeze(parsweep.betaPowBank(i,:,:)),X,Y,'spline');
V = interp2(parsweep.Rlist,parsweep.Qlist,parsweep.plot.feat,X,Y,'spline');
V(V>1e10) = NaN;
imagesc2(Xr,Xq,V)
hold on
[m,C] = contour(Xr,Xq,V,parsweep.contmap,'ShowText','on')
C.LineColor = 'k';

scatter(parsweep.InvertXY(1),parsweep.InvertXY(2),512,'kx','LineWidth',2.5)

set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
