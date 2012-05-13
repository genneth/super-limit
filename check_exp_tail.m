function check_exp_tail(pks, lambda)

n = numel(pks)-1;
norm = sum(pks .* gamma(1:(n+1)) .* lambda.^-(1:(n+1)));
mean = sum(pks .* gamma(2:(n+2)) .* lambda.^-(2:(n+2)));

fprintf('normalisation: %.2g\tmean: %.2g\n', norm, mean);

end