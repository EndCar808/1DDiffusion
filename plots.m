exact = readmatrix('exact.txt');
solution = readmatrix('solution.txt');
variables = readmatrix('variables.txt');

exact = exact(:, 1:end-1);
solution = solution(:, 1:end-1);

xmax = variables(1);
dx = variables(2);
dt = variables(3);
a = variables(4);
T = variables(5);

x = linspace(0, 2*pi, xmax);
time = 0.1:dt:T-dt;

figure
plot(x, exact(1, :), x, solution(1, :))
legend('exact', 'solution')

t = 0
figure
for i = 1:length(0.1:dt:T)-1
    t = t + dt;
    test = sin(x).*exp(-a*t);
    plot( x, solution(i, :), '-*',  x, exact(i, :), 'k')
    error(i) = norm(exact(i, :) - solution(i, :), 1);
    ylim([-1, 1])
    legend( 'exact', 'solution')
    grid on
    drawnow
end
figure
semilogy(time, error)


j = 0;
while(j < 1.0)
    
    j = j + 0.1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUDA test

cudaSolution = readmatrix('cudaSoln.txt');
cudaSolution = cudaSolution(:, 1:end-1);
figure
plot(x, cudaSolution, x , exact(1, :));
