function A = dexpinvSO3(v,input)

    %both v and input are 3x1 vectors corresponding to elements in so(3) 
    B = hat(v);
    alpha = norm(v,2);
    if alpha~=0
        dexpinvB = eye(3) - 0.5 * B + ( 1-alpha * cot(alpha/2)/2 )/(alpha^2)*B*B;
    else
        dexpinvB = eye(3);
    end
        
    A = dexpinvB * input;
end