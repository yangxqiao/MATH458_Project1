% conjugate_gradient_method
N = 3;
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
X = zeros(N*N, 1);
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

x = b;
r = b - A*x;
if norm(r) < tol
    display(x);
end
y = -r;
z = A*y;
s = y'*z;
t = (r'*y)/s;
x = x + t*y;

for k = 1:numel(b)
   r = r - t*z;
   if( norm(r) < tol )
        break;
   end
   B = (r'*z)/s;
   y = -r + B*y;
   z = A*y;
   s = y'*z;
   t = (r'*y)/s;
   x = x + t*y;
end

fprintf('The final answer obtained after %g iterations is  \n', k);
display(x)
