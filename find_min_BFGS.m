function [X_store, f_store] = find_min_BFGS(f_fun,gradf_fun,X0,X2x,x2X,n_steps)

X0 = reshape(X0,[],1);  % column vector
x0 = reshape(X2x(X0),[],1);

N = max(size(X0));
n = max(size(x0));
B0 = eye(n);

gradf0 = gradf_fun(X0);           % small
gradf0 = reshape(gradf0,[],1);    % column vector

X_store = zeros(N,n_steps+1);
X_store(:,1) = X0;

f_store = zeros(1,n_steps+1);
f_store(1) = f_fun(X0);


Xk = X0;
Bk = B0;
gradfk = gradf0;

for i = 1:n_steps

  tic
  [Xkp1, Bkp1, gradfkp1,fkp1] = BFGS_step(Xk, Bk, gradfk, f_fun,gradf_fun,x2X);  
  
  
  f_store(:,i+1) = fkp1;
  X_store(:,i+1) = Xkp1;
  
  Xk = Xkp1;
  
  Bk = Bkp1;
  gradfk = gradfkp1;
  
  tim = toc();
  fprintf('step %i done in %f s\n',i,tim)
end


end

function [Xkp1, Bkp1, gradfkp1, fkp1] = BFGS_step(Xk, Bk, gradfk, f_fun,gradf_fun,x2X)

  plot_xgradf(Xk,gradfk);


  pk = Bk\(-gradfk);
  pk = pk-linspace(pk(1),pk(end),max(size(pk)))';
  pk = pk/max(abs(pk));
  
  Pk = reshape(x2X(pk),[],1);
  Pk = Pk-linspace(Pk(1),Pk(end),max(size(Pk)))';
  Pk = Pk/max(abs(Pk));

  linsearch_options = optimoptions("fminunc",'Display','iter-detailed','OptimalityTolerance',2e-6,'MaxFunctionEvaluations',50);
  cost_fun = @(a) f_fun(Xk+a*Pk);
  [ak,~] = fminunc(cost_fun,eps,linsearch_options);

  Sk = ak*Pk;
  sk = ak*pk;
  
  Xkp1 = Xk + Sk;
  fkp1 = f_fun(Xkp1);
  
  gradfkp1 = gradf_fun(Xkp1);
  gradfkp1 = reshape(gradfkp1,[],1);    % column vector
  
  yk = gradfkp1 - gradfk;

  
  Bkp1 = Bk + yk*yk'/(yk'*sk) - Bk*(sk*sk')*Bk'/(sk'*Bk*sk);
  
 
 
end



function plot_xgradf(X,gradJ)
figure(20);
subplot(2,1,2)
plot(gradJ)
xlim([-inf inf])
hold on
subplot(2,1,1)
plot(X)
xlim([-inf inf])
hold on
drawnow
end




