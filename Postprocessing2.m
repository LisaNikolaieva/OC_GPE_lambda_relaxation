% [Psi_store] = Psi_xt(u0,grid,par,V,lambda_store(:,1));
[Psi_store] = Psi_xt(u0,grid,par,V,lambda_store(:,end));

[Psi_store2] = Psi_xt(Psi_store(:,end),grid2,par,V,zeros(1,grid2.Nt));




%%
pd = max(max(max(abs(Psi_store).^2)),max(max(abs(Psi_store2).^2)))

%%
figure
subplot(1,2,1)
imagesc([grid.t],grid.x,[abs(Psi_store).^2])
caxis([0,pd])

subplot(1,2,2)
imagesc([grid.t],grid.x,[angle(Psi_store)])
caxis([0,pd])

%%
figure
subplot(1,2,1)
imagesc([grid2.t],grid.x,[abs(Psi_store2).^2])
caxis([0,pd])

subplot(1,2,2)
imagesc([grid2.t],grid.x,[angle(Psi_store2)])
caxis([0,pd])
%%
% figure
% subplot(1,2,1)
% imagesc([grid.t,grid2.t],grid.x,[abs(Psi_store).^2,abs(Psi_store2).^2])
% caxis([0,pd])
% 
% subplot(1,2,2)
% imagesc([grid.t,grid2.t],grid.x,[angle(Psi_store),angle(Psi_store2)])
% caxis([0,pd])


%%
if 1
%% plot density
    for i = 1:50:grid.Nt
        figure(26)
        subplot(2,1,1)
        plot(grid.x,abs(Psi_store(:,i)).^2,grid.x,den_T,'--',grid.x,den_0,'--')
        grid on
        ylim([0, pd])
%         title(sprintf('linear: t = %.0f ms',grid.t(i)))
        drawnow
        subplot(2,1,2)
        plot(grid.x,angle(Psi_store(:,i)),grid.x,den_T,'--',grid.x,den_0,'--' )
        grid on
%         ylim([0, pd])
%         title(sprintf('optimized: t = %.0f ms',grid.t(i)))
        

    end
    
    
    
        for i = 1:50:grid2.Nt
        figure(26)
        subplot(2,1,1)
        plot(grid2.x,abs(Psi_store2(:,i)).^2,grid2.x,den_T,'--',grid2.x,den_0,'--')
        grid on
        ylim([0, pd])
%         title(sprintf('linear: t = %.0f ms',grid.t(i)))
        drawnow
        subplot(2,1,2)
        plot(grid2.x,angle(Psi_store2(:,i)),grid2.x,den_T,'--',grid2.x,den_0,'--' )
        grid on
%         ylim([0, pd])
%         title(sprintf('optimized: t = %.0f ms',grid.t(i)))
        

    end

end
