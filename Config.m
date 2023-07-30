close all
Physical_parameters;

grid.L = 160/2;
grid.T = 0.01;%0.001; 0.01; 0.1
T_relaxation = 1;%4.4*1;
a = 0.2; % state
b = 0.8; % energy
%% grid
grid.N=4*512;
grid.x = linspace(-grid.L/2,grid.L/2,grid.N)';
grid.dx = grid.x(2)-grid.x(1);
grid.lap4 = gallery( 'toeppen', grid.N, - 1,  16, - 30, 16, - 1 ) / ( 12 * grid.dx ^ 2 );
grid.wav = 2 * pi * ( 0 : grid.N - 1 )' / grid.N;
grid.ilap = - 2 * ( 1 - cos( grid.wav ) ) / grid.dx ^ 2;


%% define potential

Vmax        = 2e3;
sigma       = 2;
r0 = 40;
V0      =  (Vmax - Vmax * 1/2 * (1 + erf((grid.x + r0)/sigma)) + Vmax * 1/2 * ( 1+ erf((grid.x - r0)/sigma)));


%% calculate groundstates

u0 = itp(V0,par,grid); 

uT = itp_soliton(V0,par,grid);



set_V;
dV_dlambda_fun  = @(lambda) 1/(1e-3)*(V(lambda+5e-4)-V(lambda-5e-4));
mu_fun = @(lambda_t0,Psi_xt0) 0*4.3242 + 0*3.3822 + real(...
   sum(( -0.5*grid.lap4/par.m + spdiag(V(lambda_t0)+par.g*(abs(Psi_xt0).^2)) )*Psi_xt0.*conj(Psi_xt0) )...
   /sum(abs(Psi_xt0.^2)) );


mu_fun(lambda_T,uT) - mu_fun(lambda_0,u0)
%% estimate speed of sound
den_0 = abs(u0).^2;
den_T = abs(uT).^2;
figure
plot(grid.x,den_T)
hold on
plot(grid.x,0.1*max(den_T)*ones(size(den_T)),'--')
drawnow

c_s1 = sqrt(par.g*max(abs(uT).^2)/par.m); %estimated speed of sound for final state
c_s0 = sqrt(par.g*max(abs(u0).^2)/par.m); %estimated speed of sound for initial state

dt = min(0.1*grid.dx/max(c_s1,c_s0),grid.T/1000);



%% time grid
grid.Nt = ceil(grid.T/dt);
grid.t = linspace(0,grid.T,grid.Nt);
grid.dt = grid.t(2) - grid.t(1);
[grid.t_mesh,grid.x_mesh] = meshgrid(grid.t,grid.x);


%% pgrid
pgrid.L = grid.L;
pgrid.N = grid.N;
pgrid.x = grid.x;
pgrid.dx = grid.dx;
pgrid.lap4 = grid.lap4;
pgrid.nc_factor = 5;
pgrid.T = grid.T;
pgrid.Nt = ceil(grid.Nt/pgrid.nc_factor);
pgrid.t = linspace(0,pgrid.T,pgrid.Nt);
pgrid.dt = pgrid.t(2) - pgrid.t(1);         
pgrid.lap4_t   =  gallery( 'toeppen', pgrid.Nt, - 1,  16, - 30, 16, - 1 ) / ( 12 * pgrid.dt ^ 2 );
[pgrid.t_mesh,pgrid.x_mesh] = meshgrid(pgrid.t,pgrid.x);


%% grid2 & pgrid2


grid2.L = grid.L;
grid2.N=grid.N;
grid2.x = grid.x;
grid2.dx = grid.dx;
grid2.lap4 = grid.lap4;
grid2.ilap = grid.ilap;

dt = 0.1*grid2.dx/max(c_s1,c_s0);
grid2.Nt = ceil((T_relaxation-grid.T)/dt);
grid2.t = linspace(grid.T,T_relaxation,grid2.Nt);
grid2.dt = grid2.t(2) - grid2.t(1);
[grid2.t_mesh,grid2.x_mesh] = meshgrid(grid2.t,grid2.x);

pgrid2.L = grid2.L;
pgrid2.N = grid2.N;
pgrid2.x = grid2.x;
pgrid2.dx = grid2.dx;
pgrid2.lap4 = grid2.lap4;
pgrid2.nc_factor = 5;
% pgrid2.T = grid2.T;
pgrid2.Nt = ceil(grid2.Nt/pgrid2.nc_factor);
pgrid2.t = linspace(grid.T,T_relaxation,pgrid2.Nt);
pgrid2.dt = pgrid2.t(2) - pgrid2.t(1);         
pgrid2.lap4_t   =  gallery( 'toeppen', pgrid2.Nt, - 1,  16, - 30, 16, - 1 ) / ( 12 * pgrid2.dt ^ 2 );
[pgrid2.t_mesh,pgrid2.x_mesh] = meshgrid(pgrid2.t,pgrid2.x);




%% cost_fun_fun
ham       = @(lambda_t0)-0.5*grid.lap4/par.m+spdiag(V(lambda_t0));
energy_fun = @(Psi_xt0,lambda_t0) real(trapz(grid.x,conj(Psi_xt0).*(ham(lambda_t0)+par.g/2*spdiag(abs(Psi_xt0).^2))*Psi_xt0));
ET = energy_fun(uT,lambda_T);

state_fun = @(Psi_xT) 1/2*(1-abs(trapz(grid.x,conj(uT).*Psi_xT)).^2);

state_energy_fun = @(Psi_xT) a*state_fun(Psi_xT) + b*(energy_fun(Psi_xT,lambda_T)-ET);

cost = 'state_energy';
switch cost
    case 'energy'
        phi = @(Psi_xT) energy_fun(Psi_xT,lambda_T)-ET;
        pT_fun = @(Psi_xT) -2*1i*(ham(lambda_T) - spdiag(mu_fun(lambda_T,Psi_xT)*ones(size(Psi_xT))) + par.g*spdiag(abs(Psi_xT).^2))*Psi_xT;
    case 'state'
        phi = @(Psi_xT) state_fun(Psi_xT);
        pT_fun = @(Psi_xT) 1i * trapz(grid.x,conj(uT).*Psi_xT)*uT;
    case 'state_energy'
        phi = @(Psi_xT) state_energy_fun(Psi_xT);
        pT_fun = @(Psi_xT) a*(-2*1i*(ham(lambda_T) - spdiag(mu_fun(lambda_T,Psi_xT)*ones(size(Psi_xT))) + par.g*spdiag(abs(Psi_xT).^2))*Psi_xT)...
            + b*(1i * trapz(grid.x,conj(uT).*Psi_xT)*uT);
        
end

stab = @(lambda_t) par.gamma/2*sum(abs(grid.dt)*trapz((diff(lambda_t,2)/abs(grid.dt)).^2,2));
cost_fun_fun = @(Psi_xT,lambda_t) [phi(Psi_xT) + stab(lambda_t); ...
                                                    phi(Psi_xT); ...
                                                stab(lambda_t)];


int_vec_fun = @(pPsi_xt,pdVdl_xt,p_xt) trapz(pgrid.x,real(conj(pPsi_xt).*pdVdl_xt.*p_xt),1);

%% initializing
set_lambda0;

%% parameters

workspace.grid = grid;
workspace.pgrid = pgrid;
workspace.grid2 = grid2;
workspace.pgrid2 = pgrid2;
workspace.par = par;
workspace.cost_fun_fun = cost_fun_fun;
workspace.pT_fun = pT_fun;
workspace.int_vec_fun = int_vec_fun;
workspace.V = V;
workspace.dV_dlambda_fun = dV_dlambda_fun;
workspace.u0 = u0;


cf_par = workspace;
gradJ_par = workspace;
