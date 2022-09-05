function y = Forward(A, b)
    n = size(b, 1);
    y = zeros(n,1);
    y(1) = b(1) / A(1,1);
    for i = 2:n
            y(i) = (b(i) - A(i,1:(i-1))*y(1:(i-1))) / A(i,i);
    end
end