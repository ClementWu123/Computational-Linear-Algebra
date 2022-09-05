function [v,lambda,iter] = RayleighQuotient(A,v0,maxiter,tol)
    n = size(A, 1);
    v = v0;
    for k = 1:maxiter
        iter = k;
        w = (A - lambda * eye(n)) \ v;
        v = w/norm(w);
        lambda = v' * A * v;
        if norm(A * v - lambda * v) < tol
            break
        end
    end
end