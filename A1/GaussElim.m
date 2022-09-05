function x = GaussElim(A, b)
    n = size(A, 1);
    for i = 1:(n-1)
        for k = (i+1):n
            mult = A(k, i) / A(i, i);
            for j = (i+1):n
                A(k,j) = A(k,j) - mult * A(i,j);
            end
            b(k) = b(k) - mult * b(i);
        end
    end
    U = triu(A);
    x = Backward(U,b);
end