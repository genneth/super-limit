function Phat = hat(ps, z)
    n = numel(ps) - 1;
    Phat = 0;
    for k=0:n
        Phat = Phat + ps(k+1) * prod(1 + z./(1:k));
    end
end
