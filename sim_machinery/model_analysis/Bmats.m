
load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery\graph_outputs\cmap_ABC.mat')
colormap(cmap)

load('emerg_save.mat','p')
load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Rat_NPD\outputs\NPD_Final_JNPPaper_lesion\modelfit_NPD_Final_JNPPaper_lesion_2018624.mat')

R = varo;
R.chsim_name = {'MTX','STR','GPe','STN','GPI','THAL'}

p_off = R.Mfit.Pfit;

B = p_off.A{1}-p.A{1};
B(B==0) = nan;
Bi = nan(7,7);
Bi(1:end-1,1:end-1) = B;
Bmat = pcolor(Bi);
ax = gca;
ax.XTick = 1.5:1:6.5;
ax.XTickLabel = [R.chsim_name(:)]
ax.YTick = 1.5:1:6.5;
ax.YTickLabel = [R.chsim_name(:)];
 set(ax, 'YDir', 'reverse');
xlabel('From') 
ylabel('To')
c = colorbar;
c.Label.String = '\Delta P(A)';
c.Label.Rotation = -90;
c.Label.FontSize = 12;
c.Label.Position = c.Label.Position + [0.5 0 0];
title('Change in extrinsic connectivity Lesion-Control')