function p0 = calc_pT(u0,grid2,pgrid2,par,V,pT_fun)

lambda = zeros(1,grid2.Nt);
[Psi_store] = Psi_xt(u0,grid2,par,V,lambda);
plambda = interp1(grid2.t,lambda,pgrid2.t,'linear','extrap');
pPsi_store = interp2(grid2.t_mesh,grid2.x_mesh,Psi_store,pgrid2.t_mesh,pgrid2.x_mesh,"nearest",0);                       % psi_p = interp2(grid.t_mesh,grid.x_mesh,Psi_store(:,2:end),pgrid.t_mesh,pgrid.x_mesh,"nearest",0);
pT = pT_fun(Psi_store(:,end));
[p_store] = p_xt(pT,pgrid2,par,V,plambda,pPsi_store);

p0 = p_store(:,1);


        figure(21)
        subplot(2,1,1)
        plot(grid2.x,abs(Psi_store(:,end)).^2)
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');
        subplot(2,1,2)
        plot(grid2.x,angle(Psi_store(:,end)))
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');




end