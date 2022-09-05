function [A, b] = Lap2D(m)
    A = zeros(m^2);
    B = zeros(m^2, 1);
    h = 1/(m + 1);
    
    for i=1:m
        Xi = h*i;
        for j=1:m
            Yj= h * j;
            k=(j-1) * m + i;
            A(k, k) = 4;
            if (norm([Xi, Yj]-[0.35, 0.6])<=0.1 || norm([Xi, Yj]-[0.8, 0.25])<=0.1)
                B(k) = 1;
            end
            if (i > 1)
                A(k, (j - 1) * m + (i - 1)) = -1;
            end
            if (i < m)
                A(k, (j - 1) * m + (i + 1)) = -1;
            end
            if (j > 1)
                A(k, (j - 2)* m + i) = -1;
            end
            if (j < m)
                A(k, j * m + i) = -1;
            end
        end
    end
    A = A / (h^2);
    b = B;
end
