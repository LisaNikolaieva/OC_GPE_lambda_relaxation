Theta = (-1)*(grid.x<0)+(grid.x>0);

V      = @(lambda) Vmax*Theta.*lambda + V0;
