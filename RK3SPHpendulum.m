function sol = RK3SPHpendulum(vecField,action,p,h)

    k1 = zeros(6,1);
    k2 = vecField(h/2*k1,p);
    k3 = vecField(h*(-k1+2*k2), p);
    sol = action(expSE3(h*(k1/6+2*k2/3+k3/6)),p);

end