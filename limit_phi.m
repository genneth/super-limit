function [phi,phipp] = limit_phi(f, g, Tmax, umax)

m = (f(1+eps) - f(1))/eps;
alpha = find_alpha(g, m, Tmax);

if umax < 2*pi
    us = (0:0.1:umax);
elseif umax < 20*pi
    us = [(0:0.1:(2*pi)) ((2*pi):1:umax)];
else
    us = [(0:0.1:(2*pi)) ((2*pi):1:(20*pi)) ((20*pi):10:umax)];
end

% us = linspace(0, umax, 200);

% init = load('limit-phi');
% phi0 = init.phi;
phi0 = @(u) 1./(1-1.0i*u);

phi = iterate_phi(f, g, alpha, Tmax, us, phi0);

while max(abs(phi0(us) - phi(us))) > 1e-6
    fprintf('error: %.2g\n', max(abs(phi0(us) - phi(us))));
    phi0 = phi;
    phi = iterate_phi(f, g, alpha, Tmax, us, phi0);
end

fprintf('error: %.2g\n\n', max(abs(phi0(us) - phi(us))));

phi0 = phi;
[phi,phipp] = iterate_phi(f, g, alpha, Tmax, us, phi0);

% save limit-phi phi phipp

end
