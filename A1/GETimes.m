function GETimes()
    sizes = [8;16;24;32];
    gaussE = zeros(4,1);
    cholesky = zeros(4,1);
    bandGE = zeros(4,1);
    for i=1:4
        size = sizes(i);
        [A,b]=Lap2D(size);
        tic
        GaussElim(A,b);
        gaussE(i) = toc;
        tic
        Cholesky(A,b);
        cholesky(i) = toc;
        tic
        BandGE(A,b, size, size);
        bandGE(i) = toc;
    end
    t = table(sizes, gaussE, cholesky, bandGE);
    disp(t)
end