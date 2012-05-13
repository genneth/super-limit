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

end