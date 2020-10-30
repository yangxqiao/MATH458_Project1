%-------------------------------------------------------
% Runs the Jacobi method for A*x = b.
%
% ON ENTRY :
%   A        n-by-n matrix
%   b        n-dimensional vector
%   epsilon  accuracy requirement
%   maxit    maximal number of iterations
%   x        start vector for the iteration
%
% ON RETURN :
%   x        approximate solution to A*x = b.
%-------------------------------------------------------


N = 3;
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
x = zeros(N*N, 1);
epsilon = 1e-5;
maxit = 100;

for i=1:1:N
    for j=1:1:N
        A((i-1)*N+j, (i-1)*N+j) = 4;
        if i-1 > 0 && i-1 <= N
            A((i-1)*N+j, (i-2)*N+j) = -1;
        end
        if i+1 > 0 && i+1 <= N
            A((i-1)*N+j, (i)*N+j) = -1;
        end
        if j-1 > 0 && j-1 <= N
            A((i-1)*N+j, (i-1)*N+j-1) = -1;
        end
        if j+1 > 0 && j+1 <= N
            A((i-1)*N+j, (i-1)*N+j+1) = -1;
        end
    end
end

h = 1/(1+N);
for i=1:1:N
    for j=1:1:N
        if i*h <= 3/5 && i*h >= 1/5
            if j*h <= 1/2 && j*h >= 1/4
                b((i-1)*N+j, 1) = -1*h;
            end
        end
    end
end

% check if the entered matrix is a square matrix
[na, ma] = size(A);
if na ~= ma
    disp('ERROR: Matrix A must be a square matrix')
    return
end
% check if b is a column matrix
[nb, mb] = size (b);
if nb ~= na || mb~=1
   disp( 'ERROR: Matrix b must be a column matrix')
   return
end
dx = zeros(na,1);
for k=1:maxit
    sum = 0;
    for i=1:na
        dx(i) = b(i);
        for j=1:na
            dx(i) = dx(i) - A(i,j)*x(j); 
        end
        dx(i) = dx(i)/A(i,i);
        x(i) = x(i) + dx(i);
        if (dx(i) >= 0)
            sum = sum + dx(i);
        else 
            sum = sum - dx(i);
        end
    end
    if(sum <= epsilon)
        break
    end
end
fprintf('The final answer obtained after %g iterations is  \n', k);
display(x);
