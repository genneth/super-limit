function log_set

scale = arrayfun(@(n) (2*n-1) * mfun('Ei', 2*n, 2*n), 1:4);
logdist = @(n)(...
    @(x) (scale(n) * 2^(2*n-1) * n^(2*n-1) * (2*n-1)) ./ ...
         (scale(n) * x .* log(scale(n) * x).^(2*n)) ...
    );
f = @(s) s.^2;

[~,phi1pp] = limit_phi(f, logdist(1), exp(-2*1)/scale(1), 2000);
[~,phi2pp] = limit_phi(f, logdist(2), exp(-2*2)/scale(2), 2000);
[~,phi3pp] = limit_phi(f, logdist(3), exp(-2*3)/scale(3), 2000);
[~,phi4pp] = limit_phi(f, logdist(4), exp(-2*4)/scale(4), 2000);

save log-set phi1pp phi2pp phi3pp phi4pp

% phi1pp=[];phi2pp=[];phi3pp=[];%phi4pp=[];
% load log-set

umaxselector = @(t) max(min(2000.*ones(size(t)), 133./t), 133/5*ones(size(t)));
H = @(phipp, umaxs, Ws) arrayfun(@(W,umax) 1/pi * real(quadgk(@(u) fnval(phipp,u) .* exp(-1i * u * W), 0, umax)), Ws, umaxs);

t = linspace(0, 5, 500);
plot(newplot(figure) ...
    , t, H(phi1pp, umaxselector(t), t) ...
    , t, H(phi2pp, umaxselector(t), t) ...
    , t, H(phi3pp, umaxselector(t), t) ...
    , t, H(phi4pp, umaxselector(t), t) ...
    );
set(gca, 'YLim', [0 1]);

end