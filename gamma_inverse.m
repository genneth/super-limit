function gamma_inverse

inverse = @(b, k, t) (exp(-(k*b-1)*t) .* expm1(t).^((k-1)*b-1)) / beta(b, b*(k-1));

fh = figure;
ah = newplot(fh);

t = linspace(0,5,100);

set(ah, 'NextPlot', 'add');

for b=1:9
    plot(ah, t, inverse(b, 4, t), '-', 'Color', [1 0 0]);
    plot(ah, t, inverse(b, 20, t), '-', 'Color', [1 0 0]);
end

for k=2:50
    plot(ah, t, inverse(3, k, t), '-', 'Color', [0 0 1]);
end

plot(ah, t, inverse(3, 4, t), '-k', 'LineWidth', 1);
plot(ah, t, inverse(3, 20, t), '-k', 'LineWidth', 1);

set(ah, 'Box', 'off');
set(ah, 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel(ah, '$t$', 'Interpreter', 'latex', 'FontSize', 8);
ylabel(ah, '$g(t)$', 'Interpreter', 'latex', 'FontSize', 8);
 
set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', 'gamma-inverse');

close(fh);

end
