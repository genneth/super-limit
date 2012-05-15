function multi_modal

aks = [0.35 0 0  0.65];

%%
gtilde = exp_tail([0 0 1], aks); close(gcf);

umaxselector = @(t) max(min(20000.*ones(size(t)), 133./t), 133/5*ones(size(t)));
H = @(gt, umaxs, Ws) arrayfun(@(W,umax) 1/pi * real(quadgk(@(u) feval(gt,u) .* exp(-1i * u * W), 0, umax)), Ws, umaxs);

fh = figure;

ah = newplot(fh);

t = linspace(0.01, 5, 100);

powerdist = @(b)(@(t) heaviside(t-1/2) .* sinc(1/b) ./ (1+(t-1/2).^b));
pd = powerdist(2.1);
alpha = find_alpha(pd, 2, 500);

plot(ah, t, H(gtilde, umaxselector(t), t), t, pd(t/alpha)/alpha);

set(ah, 'Box', 'off');
set(ah, 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel(ah, '$t$', 'Interpreter', 'latex', 'FontSize', 8);
ylabel(ah, '$g\left(t\right)$', 'Interpreter', 'latex', 'FontSize', 8);

drawnow;

%%
% show what the distribution actually looks like
ks = 0:(numel(aks)-1);
lambda = sum(aks .* (ks+1));

ts = linspace(0,5,200);
H1 = arrayfun(@(t) sum(aks .* lambda.^(ks+1) .* t.^ks ./ gamma(ks+1)) .* exp(-lambda * t), ts);

[~,phipp] = limit_phi(@(s) s.^2, pd, 500, 20000);

inset = axes('Position', [0.5 0.5 0.4 0.4]);

plot(inset, ts, H1, ts, H(@(u) fnval(phipp, u), umaxselector(ts), ts));

set(inset, 'YLim', [0 1]);

set(inset, 'Box', 'off');
set(inset, 'FontName', 'Times New Roman', 'FontSize', 7);
xlabel(inset, '$W$', 'Interpreter', 'latex', 'FontSize', 7);
ylabel(inset, '$H\left( W \right)$', ...
     'Interpreter', 'latex', 'FontSize', 7);
 
%%
set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', 'multi-modal');

close(fh);


end