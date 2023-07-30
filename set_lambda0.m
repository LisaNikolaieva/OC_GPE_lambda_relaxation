
lambda0 = -(grid.t-grid.T/2).^2+(grid.T/2)^2;
% lambda0(1:grid.Nt/2) = grid.t(1:grid.Nt/2);
% lambda0(grid.Nt/2:grid.Nt) = -grid.t(grid.Nt/2:grid.Nt)+grid.T;


lambda0 = lambda0./max(abs(lambda0));
lambda0(1)=0;
lambda0(end)=0;