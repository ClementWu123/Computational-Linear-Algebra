function [] = findBest()
    tol = 1e-2;
    trial = 1:100;
    maxiter = 20000;
    k=8;

    size = [16, 32, 64, 128];
    alpha = [6.4e-2, 3.2e-2, 1.6e-2, 8e-3];

    for i=1:4
        m = size(i);
        a = alpha(i);
        bestOmega = 0;
        besttime = 1e8;
        bestiter=1e8;
        [~,z]=set_image(m);
        u0=FormRHS(z);
        for t=trial
            x=u0;
            totaliter=0;
            tic
            for j=1:k
                A=FormMatrix(u0,a);
                [x,iter] = SOR(t,A,u0,x,maxiter,tol);
                totaliter = totaliter+iter;
            end
            time = toc;
            if (totaliter<=bestiter && time<besttime)
                bestiter = totaliter;
                besttime = time;
                bestOmega = t;
            end
        end
        fprintf("grid size is: %d\r\n", m);
        fprintf("total time is: %d\r\n", besttime);
        fprintf("best omega is: %d\r\n\r\n", bestOmega);
    end
end