function eA = expRodrigues(w)

    A = hat(w);
    
    alpha = norm(w,2);
    
    if alpha~=0
        eA = eye(3) + sin(alpha)/alpha * A + (1-cos(alpha))/alpha^2 * A*A;
    else
        eA = eye(3) + A + 0.5 * A * A;
    end
    
end