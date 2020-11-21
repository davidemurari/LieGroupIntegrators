function A = expSE3(input)

    %we take as an input the element of SE(3) representable by (hat(u),v)
    %so u and v are column 3-vectors
    
    theta = @(w) norm(w,2);
    A = @(w) sin(theta(w))/theta(w);
    B = @(w) (1-cos(theta(w)))/(theta(w)^2);
    C = @(w) (1-A(w))/(theta(w)^2);
    
    u = input(1:3);
    v = input(4:6);
    V = @(w) eye(3) + B(w)*hat(w) + C(w) * hat(w)^2;
     
    tol = 1e-15;
    
    if theta(u)>tol
        A = [expRodrigues(v), V(v)*u];
    else
        A = [expRodrigues(v), u];
    end
   
end