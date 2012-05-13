function gtilde = exp_tail(pks, aks)

    assert(pks(1) == 0); % no extinction!
    
    Pks2 = [0 aks]; % zP(z)
    m = numel(pks)-1;
    Qks = zeros(1,numel(Pks2)*m);
    for i = 1:m
        partial = Pks2;
        for j = 2:i
            partial = conv(Pks2,partial);
        end
        for j = 1:numel(partial)
            Qks(j) = Qks(j) + pks(i+1)*partial(j);
        end
    end
    Qks = Qks(2:end); %1/z
    
    gtilde = @(ws) arrayfun(@(w) hat(aks, -1.0i * w) ./ hat(Qks, -1.0i * w), ws);
    
    % show what the distribution actually looks like
    ks = 0:(numel(aks)-1);
    lambda = sum(aks .* (ks+1));
    
    ts = linspace(0,5,200);
    H = arrayfun(@(t) sum(aks .* lambda.^(ks+1) .* t.^ks ./ gamma(ks+1)) .* exp(-lambda * t), ts);
    plot(newplot(figure), ts, H, ts, exp(-ts));

end