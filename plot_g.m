function plot_g(gtilde)

umaxselector = @(t) max(min(20000.*ones(size(t)), 133./t), 133/5*ones(size(t)));
H = @(gt, umaxs, Ws) arrayfun(@(W,umax) 1/pi * real(quadgk(@(u) feval(gt,u) .* exp(-1i * u * W), 0, umax)), Ws, umaxs);

t = linspace(0, 5, 200);
plot(newplot(figure), t, H(gtilde, umaxselector(t), t));

s = linspace(-6, 1, 200);
plot(newplot(figure), s, exp(s) .* H(gtilde, umaxselector(exp(s)), exp(s)));

end