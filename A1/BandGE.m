function x = BandGE(A,b,p,q)
    n = size(A,1);
    for k = 1:(n-1)
        for i = (k+1):min(k+p,n)
            A(i,k)=A(i,k)/A(k,k);
        end
        for i=(k+1):min(k+p,n)
            for j = (k+1):min(k+q,n)
                A(i,j)=A(i,j)-A(i,k)*A(k,j);
            end
        end
    end
    L = eye(n) + tril(A, -1);
    U = triu(A);
    y = Forward(L,b);
    x = Backward(U,y);
end