function alpha = find_alpha(g, m, tmax)

gtilde = @(s) quadgk(@(t) feval(g, t) .* exp(-s * t), 2*eps, tmax);
fprintf('gtilde(0) = %.2g\n', gtilde(0));
fprintf('gtilde(%.2g) = %.2g\n', tmax, gtilde(tmax));
if abs(gtilde(0) - 1) < 1e-5
    warning('find_alpha:normalisation', 'g is not normalised: |1-gtilde(0)| = %f', abs(gtilde(0) - 1));
end
alpha = fzero(@(s) gtilde(s) - 1/m, [0 tmax]);

end