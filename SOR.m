% successive over-relaxation
N = 3;
omega = 1;
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
x_0 = zeros(N*N, 1);
tol = 1e-5;
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

format long;
D = diag(diag(A));
L =-tril(A,-1);
U = -triu(A,1);
a = (D-omega*L);
k = 0;
for i=1:maxit
    x = a\(((1-omega)*D + omega*U)*x_0) + omega*(a\b);
    if norm(x-x_0)<tol
        break;
    end
    x_0=x;
    k = k + 1;
end
fprintf('The final answer obtained after %g iterations is  \n', k);
display(x);
