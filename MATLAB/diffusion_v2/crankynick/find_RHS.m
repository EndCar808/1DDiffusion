function RHS = find_RHS(uOld, sigma)

x_minus = circshift(uOld, [1, 0]);

x_plus = circshift(uOld, [-1, 0]);

x_centre = uOld;

RHS = sigma * x_minus + (1 - 2.0 * sigma) * x_centre + sigma * x_plus;

end