function power_set

% note that the matlab sinc(t) = sin(pi*t)/(pi*t) for some stupid reason
% powerdist = @(b)(@(t) heaviside(t-1/2) .* sinc(1/b) ./ (1+(t-1/2).^b));
% f = @(s) s.^2;

% [~,phi1pp] = limit_phi(f, powerdist(2), 500, 20000);
% [~,phi2pp] = limit_phi(f, powerdist(3), 200, 20000);
% [~,phi3pp] = limit_phi(f, powerdist(4), 50, 20000);
% [~,phi4pp] = limit_phi(f, powerdist(5), 50, 20000);

% save power-set phi1pp phi2pp phi3pp phi4pp

phi1pp=[];phi2pp=[];phi3pp=[];phi4pp=[];
load power-set

umaxselector = @(t) max(min(20000.*ones(size(t)), 133./t), 133/5*ones(size(t)));
H = @(phipp, umaxs, Ws) arrayfun(@(W,umax) 1/pi * real(quadgk(@(u) fnval(phipp,u) .* exp(-1i * u * W), 0, umax)), Ws, umaxs);

%%

fh = figure;
ah = newplot(fh);
t = linspace(0, 5, 500);
plot(ah ...
    , t, H(phi1pp, umaxselector(t), t) ...
    , t, H(phi2pp, umaxselector(t), t) ...
    , t, H(phi3pp, umaxselector(t), t) ...
    , t, H(phi4pp, umaxselector(t), t) ...
    );
set(ah, 'YLim', [0 1.5]);

set(ah, 'Box', 'off');
set(ah, 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel(ah, '$W$', 'Interpreter', 'latex', 'FontSize', 8);
ylabel(ah, '$H\left(W\right)$', 'Interpreter', 'latex', 'FontSize', 8);

legend(ah ...
    , { '$n = 2$' ...
      , '$n = 3$' ...
      , '$n = 4$' ...
      , '$n = 5$' ...
      } ...
    , 'Interpreter', 'latex', 'FontSize', 8 ...
    , 'Location', 'NorthEast' ...
    );

set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', 'power-set-W');

close(fh);

%%

fh = figure;
ah = newplot(fh);

s = linspace(-4, 2, 100);
plot(ah ...
    , s, H(phi1pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , s, H(phi2pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , s, H(phi3pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , s, H(phi4pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    );
set(ah, 'YLim', [0 1.2]);

set(ah, 'Box', 'off');
set(ah, 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel(ah, '$s=\log W$', 'Interpreter', 'latex', 'FontSize', 8);
ylabel(ah, '$J(s)$', 'Interpreter', 'latex', 'FontSize', 8);

legend(ah ...
    , { '$n = 2$' ...
      , '$n = 3$' ...
      , '$n = 4$' ...
      , '$n = 5$' ...
      } ...
    , 'Interpreter', 'latex', 'FontSize', 8 ...
    , 'Location', 'NorthEast' ...
    );

inset = axes('Position', [0.25 0.6 0.3 0.3]);
s = linspace(-6, -1, 500);
semilogy(inset ...
    , -log(-s), H(phi1pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , -log(-s), H(phi2pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , -log(-s), H(phi3pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    , -log(-s), H(phi4pp, umaxselector(exp(s)), exp(s)).*exp(s) ...
    );
set(inset, 'XLim', [min(-log(-s)) max(-log(-s))]);
set(inset, 'YLim', [10^-5 1]);
set(inset, 'YTick', 10.^(-5:0));

set(inset, 'Box', 'off');
set(inset, 'FontName', 'Times New Roman', 'FontSize', 7);
xlabel(inset, '$-\log |s|$', 'Interpreter', 'latex', 'FontSize', 7);
ylabel(inset, '$J(s)$', ...
     'Interpreter', 'latex', 'FontSize', 7);

set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', 'power-set-log-W');

close(fh);

%%

end