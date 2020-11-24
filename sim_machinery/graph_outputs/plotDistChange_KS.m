function plotDistChange_KS(Rho,nu,xf,psave,pInd,R,stdev)
% Univariate plots
indFlat = spm_vec(pInd);
subplot(2,2,1)
for i = 1:2
    if i == 1
        ls = '--';
        nind = 1;
        
        %{'.params','.noisecov'}
        %         M = psave(nind).params; M(M==0) = [];
        if strcmp(R.projectn,'MVAR')
            M = psave(nind).params; M_s = psave(nind).params_s;
            M_s(M(:)==0) = []; M(M(:)==0) = [];
            Ma_s = M_s(M(:)>-30); Ma = M(M(:)>-30);
            
        elseif strcmp(R.projectn,'Intracranial_pressure')
            M = [psave(nind).Ca psave(nind).Ro psave(nind).Can psave(nind).KE psave(nind).Pc];
            M_s = [psave(nind).Ca_s psave(nind).Ro_s psave(nind).Can_s psave(nind).KE_s psave(nind).Pc_s];
            %M_s(M(:)==0) = [];% M(M(:)==0) = [];
            Ma_s = M_s(M(:)>-30); Ma = M(M(:)>-30);
        else
            M = psave(nind).A{1}; M_s = psave(nind).A_s{1};
            %M_s(M(:)==0) = [];% M(M(:)==0) = [];
            Ma_s = M_s(M(:)>-30); Ma = M(M(:)>-30);
        end
        cmap = linspecer(5);
        X = R.SimAn.pOptBound(1):.05:R.SimAn.pOptBound(2);
        QL = length(Ma);
        if QL>5
            QL = 5;
        end
        for Q = 1:QL
            p = normpdf(X,Ma(Q),Ma_s(Q).*R.SimAn.jitter);
            %             [p,type,coefs] = pearspdf(X,Ma(Q),R.SimAn.jitter*R.SimAn.Tm,1,3);
            plot(X,p,ls,'color',cmap(Q,:))
            hold on
        end
    else % Copula
        ls = '-';
        r = copularnd('t',Rho,nu,500);
        for Q = 1:5 %size(xf,1)
            x1 = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
            [x1,f] = ksdensity(x1,R.SimAn.pOptRange);
            plot(f,x1,ls,'color',cmap(Q,:))
            hold on
        end
    end
end
xlabel('\mu')
ylabel('p(\mu)')
ylim([0 1]);
xlim(R.SimAn.pOptBound.*0.2)
title('Approximate Posterior Distribution')
axis square

subplot(2,2,2)
imagesc(Rho)
set(gca,'YDir','normal')
title('Copula Covariance')
set(gca,'XTick',1:size(Rho,1))
set(gca,'YTick',1:size(Rho,1))
colormap(brewermap(128,'RdYlBu'))


% Multivariate plots
axis square
r = copularnd('t',Rho,nu,1000);
subplot(2,2,3)
title('2D Sample Drawn from Copula')
if strcmp(R.projectn,'MVAR')
    i = pInd.params(1);
    j =pInd.params(2);
else
    i = spm_vec(pInd);
    i = i(4);
    j = spm_vec(pInd);
    j = j(5);
end
i = find(indFlat==i);
j = find(indFlat==j);
u1 = r(:,i);
v1 = r(:,j);
x1 = ksdensity(xf(i,:),u1,'function','icdf');
y1 = ksdensity(xf(j,:),v1,'function','icdf');
scatter(x1,y1);
xlim(R.SimAn.pOptBound.*0.2)
ylim(R.SimAn.pOptBound.*0.2)

xlabel('M1 Time Constant'); ylabel('M1 Synaptic Gain')
set(get(gca,'children'),'marker','.')
axis square
subplot(2,2,4)
title(' 3D Sample Drawn from Copula')
if strcmp(R.projectn,'MVAR')
    i = pInd.params(1);
    j =pInd.params(2);
    k = pInd.params(3);
else
    i = spm_vec(pInd);
    i = i(4);
    j = spm_vec(pInd);
    j = j(5);
    k = spm_vec(pInd);
    k = k(2);
end
i = find(indFlat==i);
j = find(indFlat==j);
k = find(indFlat==k);

u1 = r(:,i);
v1 = r(:,j);
w1 = r(:,k);
x1 = ksdensity(xf(i,:),u1,'function','icdf');
y1 = ksdensity(xf(j,:),v1,'function','icdf');
z1 = ksdensity(xf(k,:),w1,'function','icdf');

scatter3(x1,y1,z1)
xlim(R.SimAn.pOptBound.*0.2)
ylim(R.SimAn.pOptBound.*0.2)
zlim(R.SimAn.pOptBound.*0.2)

xlabel('M1->STR A'); ylabel('STR Gain'); zlabel('M1->STR D');
set(get(gca,'children'),'marker','.')
axis square
set(gcf,'Position',[2.5000  272.0000  903.0000  728.0000])
