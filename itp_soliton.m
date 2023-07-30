function Psi = itp_soliton(V,par,grid)

mue_0 = 0;
n = 3000;


v = 0;
gamma = 1/sqrt(1 - v^2);
xi_s = 1;%4*grid.dx;
cs = 1;
t=1;
f_tanh = (grid.x - (v*cs*t))/(sqrt(2)*gamma*xi_s);

Psi = (1j*v + (1/gamma)*tanh(f_tanh));%*exp(1j*mu*t);





normalize = @(Psi)Psi/sqrt(trapz(grid.x,conj(Psi).*Psi));
% Psi = normalize(Psi);
for i = 1:n
        Psi_next = propagate_Psi_SS(-1i*0.01,grid,par,Psi,V,mue_0);
    Psi_next = normalize(Psi_next);
    
    if trapz(grid.x,abs(Psi_next - Psi))<1e-9
        break;
    end
    Psi = Psi_next;
end


end

