function I = fct_integrand_noise_energy(model,kx,ky)
% Compute the integrand involved in the computation of M
%

alpha = (3-model.sigma.slope_sigma)/2;

I(1,1,:,:) = 1/(2*(2-alpha)) * kx.^(3-2*alpha) .* ky .* ...
    inbeta( - (ky ./ kx).^2 , 1/2 , 2-alpha);
%     fct_beta_inc( - (ky ./ kx).^2 , 1/2 , 2-alpha);

k2 = kx.^2+ky.^2;
I(1,2,:,:) = 1/(4*(1-alpha)*(2-alpha)) * k2.^(2-alpha);

end

% function bet = fct_beta_inc(x,p,q)
% % Incomplete Beta function without normalization
% %
% 
% bet = beta(p,q) * betainc(x,p,q);
% end