function f_pos = fct_smooth_pos_part(f)
% Compute a smooth approximation of the positive part fo f (i.e. the
% function which is equal to f where f > 0 and equal to zero where f is
% negative)
%

sf = size(f); f=f(:);

sigma_f = sqrt(mean(f(:).^2));
order = 8.;
coef_sigma = 3;
sigma_f = sigma_f / coef_sigma;
f_pos = (1 - exp(-(1/sigma_f * f ).^order)) ;
f_pos( f < 0) = 0  ;
f_pos = f_pos .* f;

f_pos = reshape(f_pos,sf);