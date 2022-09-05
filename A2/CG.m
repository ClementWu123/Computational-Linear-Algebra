function [x,iter] = CG(A,b,x_initial,maxiter,tol)
    x = x_initial;
    r = b - A*x;
    p = r;
    rnew = r'*r;
    for iter = 0:maxiter
        if (norm(r)<=tol*norm(b))
            break
        end
        alpha = rnew / (p'*A*p);
        x = x + alpha * p;
        r = r - alpha * A * p;
        rold = rnew;
        rnew = r' * r;
        beta = rnew / rold;
        p = r + beta * p;
    end
end