function gtilde = exp_tail(aks, pks, lambda)

    assert(aks(1) == 0); % no extinction!
    
    ks = 0:(numel(pks) - 1);
    Pks = gamma(1+ks) .* pks;
    
    Pks2 = [0 Pks]; % zP(z)
    m = numel(aks)-1;
    Qks = zeros(1,numel(Pks2)*m);
    for i = 1:m
        partial = Pks2;
        for j = 2:i
            partial = conv(Pks2,partial);
        end
        for j = 1:numel(partial)
            Qks(j) = Qks(j) + aks(i+1)*partial(j);
        end
    end
    Qks = Qks(2:end); %1/z
    
    gtilde = @(ws) arrayfun(@(w) hat(Pks, lambda, -1.0i * w) ./ hat(Qks, lambda, -1.0i * w), ws);

end