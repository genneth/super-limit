function [phi1,phi1pp] = iterate_phi(f, g, alpha, Tmax, us, phi0)

assert(us(1) == 0);

ps = zeros(size(us));

ps(1) = 1;

for i = 2:numel(us) % notice that we skip zero

    ps(i) = quadl(@(y) feval(g,y) .* feval(f, feval(phi0, us(i) * exp(-alpha*y))), 0, Tmax, 1e-8);

end

phi1pp = csape(us, [1i ps 0], [1 2]); % fix first derivative at 0, and second derivative is zero at far end
phi1 = @(x) fnval(phi1pp, x);

end
