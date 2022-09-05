function [] = denoise()
    tol=1e-2;
    maxiter=20000;
    k=8;
    size = [16,32,64,128];
    alpha = [6.4e-2, 3.2e-2, 1.6e-2, 8e-3];
    omega = [14, 87, 79, 88];
    image_size = zeros(16,1);
    method = strings(16,1);
    times = zeros(16,1);
    iterations = zeros(16,1);

    for i=1:4
        m = size(i);
        a = alpha(i);
        o = omega(i);
        [~,z]=set_image(m);
        u0 = FormRHS(z);
        iterjacobi = 0;
        itergs = 0;
        itersor = 0;
        itercg = 0;
        x_j = u0;
        x_g = u0;
        x_s = u0;
        x_c = u0;

        method((i-1)*4+1)="jacobi";
        image_size((i-1)*4+1)=m;
        tic
        for j=1:k
            A = FormMatrix(x_j,a);
            [x_j,iter] = Jacobi(A,u0,x_j,maxiter,tol);
            iterjacobi = iterjacobi + iter;
        end
        time = toc;
        times((i-1)*4+1)=time;
        iterations((i-1)*4+1)=iterjacobi;

        method((i-1)*4+2)="GS";
        image_size((i-1)*4+2)=m;
        tic
        for j=1:k
            A = FormMatrix(x_g,a);
            [x_g,iter] = GS(A,u0,x_g,maxiter,tol);
            itergs = itergs + iter;
        end
        time = toc;
        times((i-1)*4+2)=time;
        iterations((i-1)*4+2)=itergs;

        method((i-1)*4+3)="SOR";
        image_size((i-1)*4+3)=m;
        tic
        for j=1:k
            A = FormMatrix(x_s,a);
            [x_s,iter] = SOR(o,A,u0,x_s,maxiter,tol);
            itersor = itersor + iter;
        end
        time = toc;
        times((i-1)*4+3)=time;
        iterations((i-1)*4+3)=itersor;

        method((i-1)*4+4)="CG";
        image_size((i-1)*4+4)=m;
        tic
        for j=1:k
            A = FormMatrix(x_c,a);
            [x_c,iter] = CG(A,u0,x_c,maxiter,tol);
            itercg = itercg + iter;
        end
        time = toc;
        times((i-1)*4+4)=time;
        iterations((i-1)*4+4)=itercg;
    end
    t = table(method,image_size,times,iterations);
    disp(t);

end