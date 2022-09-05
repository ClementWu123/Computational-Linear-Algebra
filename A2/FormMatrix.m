function A = FormMatrix(u, alpha)

    beta = 1e-6;
    n = size(u, 1);
    m = sqrt(n);
    v = reshape(u,[m,m])';
    h = 1/(1+m); 

    aw = zeros(m);
    ae = zeros(m);
    as = zeros(m);
    an = zeros(m);
    ac = zeros(m);

    function s = get(i ,j)
        if (i>m || i<1 || j>m || j<1)
            s = 0;
        else
            s = v(i ,j);
        end
    end

    function s = calc(v1,v2,v3,v4)
        s = 1/(2*sqrt(((v1-v2)/h)^2+((v3-v4)/h)^2+beta));
    end

    for i = 1:m
        for j = 1:m
            aw(i,j) = (-alpha/h^2)*(calc(get(i,j), get(i-1,j), get(i,j), get(i,j-1))+...
                calc(get(i,j), get(i-1,j), get(i-1,j+1), get(i-1,j)));
            ae(i,j) = (-alpha/h^2)*(calc(get(i+1,j), get(i,j), get(i+1,j), get(i+1,j-1))+...
                calc(get(i+1,j), get(i,j), get(i,j+1), get(i,j)));
            as(i,j) = (-alpha/h^2)*(calc(get(i,j), get(i-1,j), get(i,j), get(i,j-1))+...
                calc(get(i+1,j-1), get(i,j-1), get(i,j), get(i,j-1)));
            an(i,j) = (-alpha/h^2)*(calc(get(i+1,j), get(i,j), get(i,j+1), get(i,j))+...
                calc(get(i,j+1), get(i-1,j+1), get(i,j+1), get(i,j)));
            ac(i,j) = -(aw(i,j)+ae(i,j)+as(i,j)+an(i,j))+1;
        end
    end

    wcol = reshape(aw(2:m,:)',[],1);
    ecol = reshape(ae(1:m-1,:)',[],1);
    scol = reshape(as',[],1);
    ncol = reshape(an',[],1);
    ccol = reshape(ac',[],1);
    A = sparse(diag(wcol,-m)+diag(ecol,m)+diag(scol(2:n),-1)+diag(ncol(1:n-1),1)+diag(ccol));
end