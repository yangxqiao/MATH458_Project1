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


N = 30;
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
x = zeros(N*N, 1);
epsilon = 1e-5;
maxit = 10000;

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

Norm = [];
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
    
    temp = norm(b-A*x, 2);
    Norm = [Norm, temp];
end
fprintf('The final answer obtained after %g iterations is  \n', k);
display(x);
plot(Norm, 'r');

hold on;


% conjugate_gradient_method
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
X = zeros(N*N, 1);
tol = 1e-5;
maxit = 10000;

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

Norm = [];
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
   
   temp = norm(b-A*x, 2);
   Norm = [Norm, temp];
end

fprintf('The final answer obtained after %g iterations is  \n', k);
display(x)
plot(Norm, 'b');

hold on;

% successive over-relaxation
omega = 1.172;
A = zeros(N*N, N*N);
b = zeros(N*N, 1);
x_0 = zeros(N*N, 1);
tol = 1e-5;
maxit = 10000;

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

Norm = [];
for i=1:maxit
    x = a\(((1-omega)*D + omega*U)*x_0) + omega*(a\b);
    if norm(x-x_0)<tol
        break;
    end
    x_0=x;
    k = k + 1;
    
    temp = norm(b-A*x, 2);
    Norm = [Norm, temp];
end
fprintf('The final answer obtained after %g iterations is  \n', k);
display(x);
plot(Norm, 'g');
