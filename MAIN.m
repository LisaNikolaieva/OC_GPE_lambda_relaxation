close all;
Config;

%%

f_fun = @(X) cost_function(X,cf_par); % f:      X |-> R
gradf_fun = @(X) gradJ(X,gradJ_par);  % gradf:  X |-> x
X0 = lambda0;
X2x = @(X) interp1(grid.t,X,pgrid.t,'linear','extrap'); % X2x:  X |-> x
x2X = @(x) interp1(pgrid.t,x,grid.t,'linear','extrap'); % x2X:  x |-> X
n_steps = 10;

[lambda_store, cost_function_store] = find_min_BFGS(f_fun,gradf_fun,X0,X2x,x2X,n_steps);
%%
Postprocessing2;
