function A = hat(v)

    A = zeros(3);
    A(1,2) = -v(3);
    A(1,3) = v(2);
    A(2,3) = -v(1);
    
    A = A - A';
    
end