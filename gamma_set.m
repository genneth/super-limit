function gamma_set

%gammadist = @(b)(@(t)b^b/gamma(b)*t.^(b-1).*exp(-b*t));
% f = @(s) s.^2;
% 
% [~,phi1pp] = limit_phi(f, gammadist(1), 20, 2000);
% [~,phi2pp] = limit_phi(f, gammadist(2), 20, 2000);
% [~,phi3pp] = limit_phi(f, gammadist(3), 20, 2000);
% [~,phi4pp] = limit_phi(f, gammadist(4), 20, 2000);
% 
% save gamma-set phi1pp phi2pp phi3pp phi4pp

phi1pp=[];phi2pp=[];phi3pp=[];phi4pp=[];
load gamma-set

umaxselector = @(t) max(min(2000.*ones(size(t)), 133./t), 133/5*ones(size(t)));
H = @(phipp, umaxs, Ws) arrayfun(@(W,umax) 1/pi * real(quadgk(@(u) fnval(phipp,u) .* exp(-1i * u * W), 0, umax)), Ws, umaxs);

fh = figure;

ah = newplot(fh);

t = linspace(0.01, 5, 100);
plot(ah ...
    , t, H(phi1pp, umaxselector(t), t) ...
    , t, H(phi2pp, umaxselector(t), t) ...
    , t, H(phi3pp, umaxselector(t), t) ...
    , t, H(phi4pp, umaxselector(t), t) ...
    );
set(ah, 'YLim', [0 1]);

set(ah, 'Box', 'off');
set(ah, 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel(ah, '$W$', 'Interpreter', 'latex', 'FontSize', 8);
ylabel(ah, '$H\left(W\right)$', 'Interpreter', 'latex', 'FontSize', 8);

guess = @(b,t) b^b/gamma(b)*t.^(b-1).*exp(-b*t);

inset = axes('Position', [0.5 0.5 0.4 0.4]);
plot(inset ...
    , t, H(phi1pp, umaxselector(t), t) - guess(1,t) ...
    , t, H(phi2pp, umaxselector(t), t) - guess(2,t) ...
    , t, H(phi3pp, umaxselector(t), t) - guess(3,t) ...
    , t, H(phi4pp, umaxselector(t), t) - guess(4,t) ...
    );
set(inset, 'XLim', [0 5]);
set(inset, 'YLim', [-0.03 0.03]);

set(inset, 'Box', 'off');
set(inset, 'FontName', 'Times New Roman', 'FontSize', 7);
xlabel(inset, '$W$', 'Interpreter', 'latex', 'FontSize', 7);
ylabel(inset, 'residual from $\Gamma$-distribution', ...
     'Interpreter', 'latex', 'FontSize', 7);

legend(inset ...
    , { '$\beta = 1$' ...
      , '$\beta = 2$' ...
      , '$\beta = 3$' ...
      , '$\beta = 4$' ...
      } ...
    , 'Interpreter', 'latex', 'FontSize', 7 ...
    , 'Location', 'NorthEast' ...
    );

set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', 'gamma-set');

close(fh);

end