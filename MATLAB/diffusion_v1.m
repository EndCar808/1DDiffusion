tsteps = 60;
xsteps = 60;
XMAX   = xsteps + 1;

T = 1.0; 
t0 = 0.1;
dt = (T - t0)/tsteps;
time = t0:dt:T;

dx = 2*pi/(XMAX);
x = 0:dx:(2*pi - dx);

a = 0.1;
r = (a*dt)/(2.0*(dx^2));

b = zeros(XMAX, 1);
dd = (1 + 2*r)*ones(XMAX, 1);
dl = -r*ones(XMAX - 1, 1);
A = diag(dd, 0) + diag(dl, 1) + diag(dl, -1);


u      = ones(XMAX, 1);
v      = ones(XMAX, 1);
u(1)   = -(1 + 2*r);
u(end) = -r;
v(1)   = 1;
v(end) = -(-r/(1 + 2*r));

A(1,1)      = A(1, 1) - u(1);
A(end, end) = A(end, end) - (u(end)*v(end));



U = sin(x);
wrap = @(indx, n) (1 + mod(indx - 1, n)); % use to wrap the indexing around


figure
for t = 1:length(time)
    b(1) = r*U(end) + (1 - 2*r)*U(1) + r*U(2);
    for i = 2:XMAX-1
        b(i) = r*U(i - 1) + (1-2*r)*U(i) + r*U(i + 1);
    end
    b(end) = r*U(end - 1) + (1 - 2*r)*U(end) + r*U(1);
     
    for i = 1:XMAX
        b(i) = r*U(wrap(i - 1,XMAX)) + (1-2*r)*U(i) + r*U(wrap(i + 1, XMAX));
    end

    y = A\b;
    z = A\u;
    
    
    U = y - ((v.*y)/(1 + v.*z))*z;
    plot( x, sin(x)*exp(-a*time(t)), x, U)
    ylim([-1, 1])
    legend( 'exact', 'solution')
    grid on
    drawnow
end