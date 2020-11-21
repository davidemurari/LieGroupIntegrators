function vec = dexpinvSE3 (sigma, input)

    A = sigma(1:3);
    a = sigma(4:6);
    
    rho = A'*a;
    g1 = @(z) -0.5;
    g1tilde = @(z) 0;
    g2 = @(z) (1-z/2*cot(z/2))/(z^2);
    g2prime = @(z) 1/z^2 * (-0.5*cot(z/2)-z/2*(-0.5-0.5*cot(z/2)^2))-...
        2/z^3 * (1-z/2 * cot(z/2));
    g2tilde = @(z) rho/z * g2prime(z);
    f0 = 1;
    alpha = norm(A,2);
    mat = zeros(3);
    tol = 1e-15;
    
    if alpha>tol
        Mat = f0*eye(6) + [g1(alpha)*hat(A), mat; g1(alpha)*hat(a)+...
            g1tilde(alpha)*hat(A), g1(alpha)*hat(A)] + ... 
            [hat(A), mat;hat(a), hat(A)] * [g2(alpha)*hat(A), mat; ...
            g2(alpha)*hat(a)+g2tilde(alpha)*hat(A), g2(alpha)*hat(A)];
        vec = Mat * input;
    else
        vec = input;
    end
end