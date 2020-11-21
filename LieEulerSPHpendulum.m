function sol = LieEulerSPHpendulum(vecField,action,p,h)

    k0 = zeros(6,1);
    k1 = vecField(k0,p);
    sol = action(expSE3(h*k1),p);

end