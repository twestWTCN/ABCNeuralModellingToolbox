        


r = copularnd('t',Rho,nu,500);
        for Q = size(xf,1)
            indFlat
            
            
            x1 = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
            [x1,f] = ksdensity(x1,R.SimAn.pOptRange);
            plot(f,x1,ls,'color',cmap(Q,:))
            hold on
        end