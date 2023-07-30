function [X_store, f_store] = find_min_BFGS2(f_fun,gradf_fun,X0,X2x,x2X,n_steps)

X0 = reshape(X0,[],1);      % column vector
x0 = reshape(X2x(X0),[],1);

N = max(size(X0));
n = max(size(x0));
B0 = eye(n);

X_store = zeros(N,n_steps+1);
X_store(:,1) = X0;

f_store = zeros(1,n_steps+1);
f_store(1) = f_fun(X0);


Xkm1 = X0;
Bkm1 = B0;
skm1 = zeros(n,1);
is_initial = 1;
for i = 1:n_steps

  tic
 [Bk,sk,Xk,fk] = make_BFGS_step(Bkm1,skm1,Xkm1,f_fun,gradf_fun,x2X,is_initial);
  
 f_store(:,i+1) = fk;
 X_store(:,i+1) = Xk;
  
 Bkm1 = Bk;
 skm1 = sk; 
 Xkm1 = Xk;

  is_initial = 0;
  tim = toc();
  fprintf('step %i done in %f s\n',i,tim)
end


end
   
function [Bk,sk,Xk,fk] = make_BFGS_step(Bkm1,skm1,Xkm1,f_fun,gradf_fun,x2X,is_initial)   

n = max(size(skm1));



grad_J = reshape(gradf_fun(Xkm1), [], 1);  % column vector
plot_xgradf(Xkm1,grad_J);

ykm1 = grad_J;



if is_initial
    Bk = repmat(eye(n),[1,1,1]);
else
skm1 = reshape(skm1,[],1);
Sk = skm1*skm1';
SYk = skm1*ykm1';
yksk = ykm1'*skm1;
Bk = (speye(n)-sparse(SYk)/yksk)*Bkm1*(speye(n)-sparse(SYk')/yksk)+Sk/yksk;
end



pk = -Bk\grad_J;

pk = pk./max(abs(pk));
Pk = reshape(x2X(pk),[],1);
Pk = Pk - linspace(Pk(1),Pk(end),max(size(Pk)))';
pk = pk - linspace(pk(1),pk(end),max(size(pk)))';


cost_fun = @(a) f_fun(Xkm1 + a.*Pk);


linsearch_options = optimoptions("fminunc",'Display','iter-detailed');
[a_opt,fk] = fminunc(cost_fun,eps*ones(1,1),linsearch_options);
fprintf('optimal a: %g\n',a_opt')

sk = a_opt.*pk;

Xk = Xkm1 + a_opt.*Pk;
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




