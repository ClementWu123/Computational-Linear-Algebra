function NL = CreateImageGraph(U)

[rows, cols]=size(U);
n = rows*cols;

sigdist = 150;
sigint = 0.002;
W = zeros(n);

for i=1:rows
    for j=1:cols
       current = U(i,j);
       if (i-1)>0
           top = U(i-1,j);
           W(j+(i-1)*cols,j+(i-2)*cols) = exp(-1/sigdist)...
               *exp(-(current-top)^2/sigint);
       end
       if (i-1)>0 && (j-1)>0
           topleft = U(i-1,j-1);
           W(j+(i-1)*cols,(j-1)+(i-2)*cols) = exp(-2/sigdist)...
               *exp(-(current-topleft)^2/sigint);
       end
       if (i-1)>0 && (j+1)<=cols
           topright = U(i-1,j+1);
           W(j+(i-1)*cols,(j+1)+(i-2)*cols) = exp(-2/sigdist)...
               *exp(-(current-topright)^2/sigint);
       end
       if (j-1)>0
           left = U(i,j-1);
           W(j+(i-1)*cols,(j-1)+(i-1)*cols) = exp(-1/sigdist)...
               *exp(-(current-left)^2/sigint);
       end
       if (j+1)<=cols
           right = U(i,j+1);
           W(j+(i-1)*cols,(j+1)+(i-1)*cols) = exp(-1/sigdist)...
               *exp(-(current-right)^2/sigint);
       end
       if (i+1)<=rows 
           bottom = U(i+1,j);
           W(j+(i-1)*cols,j+i*cols) = exp(-1/sigdist)...
               *exp(-(current-bottom)^2/sigint);
       end
       if (i+1)<=rows && (j-1)>0
           bottomleft = U(i+1,j-1);
           W(j+(i-1)*cols,(j-1)+i*cols) = exp(-2/sigdist)...
               *exp(-(current-bottomleft)^2/sigint);
       end
       if (i+1)<=rows && (j+1)<=cols
           bottomright = U(i+1,j+1);
           W(j+(i-1)*cols,(j+1)+i*cols) = exp(-2/sigdist)...
               *exp(-(current-bottomright)^2/sigint);
       end
    end
end

D = diag(sum(W));
NL = sparse(eye(n) - D^(-1/2)*W*D^(-1/2));

