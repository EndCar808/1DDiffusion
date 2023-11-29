n = 128;
L = 2 * pi;
dx = L / n;
dt = 1e-03;

T = 1.0;

sigma = 0.5 * dt / (dx ^ 2);

for i = 1:n
    x(i) = (i - 1) * dx;
    uOld(i) = cos((i - 1) * dx);
end

RHS = find_RHS(uOld, sigma);

A = zeros(n, n);

a = - sigma;
b = 1.0 + 2.0 * sigma;
c = - sigma;

for i = 1:n
    A(i, i) = b;
    
    if i > 1
        A(i, i - 1) = a;
    end
    
    if i < n
        A(i, i + 1) = c;
    end
    
end

alpha = c;
beta = a;
gamma = - b;

u = zeros(n, 1);
u(1) = gamma;
u(n) = alpha;

v = zeros(n, 1);
v(1) = 1.0;
v(n) = beta / gamma;

A(1, 1) = A(1, 1) - gamma;
A(n, n) = A(n, n) - alpha * beta / gamma;

q = A \ u;

time = 0;

while time < T
    RHS = find_RHS(uOld, sigma);

    y = A \ RHS;
    
    uOld = y - ((transpose(v) * y) / (1 + transpose(v) * q)) * q;
    
    time = time + dt;
end