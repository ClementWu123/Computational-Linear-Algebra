function [V,Lambda,iter] = QRIteration(A,maxiter,tol)
    intermediate = A;
    n = size(A,1);
    V = eye(n);
    for iter = 1:maxiter
        [Q,R] = qr(intermediate);
        V = V * Q;
        intermediate = R * Q;
        Lambda = diag(intermediate);
        if all(norm(A .* V - Lambda' .* V)) < tol
            break
        end
    end
end