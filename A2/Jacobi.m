function [x,iter] = Jacobi(A,b,x_initial,maxiter,tol)
    D = diag(diag(A));
    x = x_initial;
    for iter = 0:maxiter
        r = b - A*x;
        if (norm(r)<=tol*norm(b))
            break
        end
        x = x + D^(-1)*r;
    end
end