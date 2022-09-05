function EigenMethods()
    n = 100;
    A = diag(2*ones(1,n))+diag(-1*ones(1,n-1),1)+diag(-1*ones(1,n-1),-1);
    maxiter = 10000;
    tol = 1e-4;
    vp0 = zeros(n,1);
    vp0(1) = 1;
    vr0 = ones(n,1);

    [vp,lambdap,iterp] =  PowerIteration(A,vp0,maxiter,tol);
    plot(vp);
    title(...
        sprintf("Power Interation - Largest Eigenvalue %f at %d iteration",...
        lambdap, iterp));
    saveas(gcf,'PowerIteration.png');

    [vr,lambdar,iterr] =  PowerIteration(A,vr0,maxiter,tol);
    plot(vr);
    title(...
        sprintf("Rayleigh Quotient - Largest Eigenvalue %f at %d iteration",...
        lambdar, iterr));
    saveas(gcf,'RayleighQuotient.png');

    [V,Lambda,iter] = QRIteration(A,maxiter,tol);
    plot(Lambda);
    title("All Eigenvalues of QRIteration");
    saveas(gcf,'AllEigenvalues.png');

    for i = 1:4
        pos = 20 * i;
        plot(V(:,pos));
        title(...
        sprintf("QRIteration - Number %d Eigenvector with ..." + ...
        "Eigenvalue %f at %d iter", pos, Lambda(pos), iter));
        saveas(gcf,sprintf("QRIteration%d.png",pos));
    end

end