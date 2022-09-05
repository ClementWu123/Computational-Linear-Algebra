function [x,iter] = GS(A,b,x_initial,maxiter,tol)
    D = diag(diag(A));
    L = tril(A);
    x = x_initial;
    for iter = 0:maxiter
        r = b - A*x;
        if (norm(r)<=tol*norm(b))
            break
        end
        x = x + (D+L) \ r;
    end
end