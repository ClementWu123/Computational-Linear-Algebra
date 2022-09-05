function [v,lambda,iter] =  PowerIteration(A,v0,maxiter,tol)
    v = v0;
    for k = 1:maxiter
        iter = k;
        w = A*v;
        v = w/norm(w);
        lambda = v' * A * v;
        if norm(A * v - lambda * v) < tol
            break
        end
    end
end
