% Question 3: experiment with relaxation parameter for SOR
% construct matrices A, U, D, and L for SOR method
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

U = zeros(9);
D = zeros(9);
L = zeros(9);
for i = 1:9
    for j = 1:9
        if i < j
            L(i,j) = A(i,j);
        elseif i == j
            D(i,j) = A(i,j);
        else
            U(i,j) = A(i,j);
        end
    end
end

% experiment with w values
row = 1;
eigVec = zeros(1,2001);
for w=0:0.001:2
    intermediate1 = inv(L.*w+D);
    intermediate2 = D.*(1-w) - w.*U;
    T = intermediate1 * intermediate2;
    e = eig(T);
    emax = abs(max(e));
    eigVec(1,row) = emax;
    row = row + 1;
end
wVec = 0:0.001:2;
plot(wVec,eigVec,'-r')



% function for calculating g(x,y)
function g = g(x,y)
    if(x>=0.2 && x<= 0.6 && y>= 0.25 && y<=0.5)
        g = 1;
    else
        g = 0;
    end
end

