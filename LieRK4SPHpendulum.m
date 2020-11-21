function sol = LieRK4SPHpendulum(action,F,p,h)
    %Runge kutta Heun
   
    field = @(sigma,p) action(expSE3(sigma),p);
    
    k1 = h*F(p);
    Y2 = field(k1/2,p);
    k2 = h*F(Y2);
    k3 = h*F(field(k2/2,p));
    k4 = h*F(field(k3-k1/2,Y2));
    yHalf = field(1/12 * (3*k1+2*k2+2*k3-k4),p);
    
    
    sol = field(1/12 * (-k1+2*k2+2*k3+3*k4),yHalf);
    
end