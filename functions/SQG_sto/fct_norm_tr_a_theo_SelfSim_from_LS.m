function norm_tr_a_theo = fct_norm_tr_a_theo_SelfSim_from_LS(...
    model,kappa1,kappa2)
% Compute the theoretical trace of the variance tensor a
%
error('This function should not be used');

alpha = ( 3 - model.sigma.slope_sigma )/2;

% %d_kappa = (2*pi) / sqrt(prod(model.grid.MX.*model.grid.dX));
% dkxdky = (2*pi)^2 / (prod(model.grid.MX.*model.grid.dX));

% Calcul of energy
% One has to divid by prod(model.grid.MX)^2 because of the 2 inverse
% Fourier transform.
% Moreover, one has to divide by d_kappa in order to transform the factor
% dk_x d_ky (necessary for the analitic integration) into d_kappa which
% appear in the definition of sigma_on_sq_dt

% % % %norm_a_theo = (pi/(4-2*alpha)) * ...
% % % %norm_a_theo = (1/(4-2*alpha)) * ...
% % % norm_tr_a_theo = (1/( prod(model.grid.MX) * d_kappa )) * ...
% % norm_tr_a_theo = (1/( prod(model.grid.MX)^2 * dkxdky )) * ...
% norm_tr_a_theo = (1/ dkxdky ) * ...
%     (2*pi/(4-2*alpha)) * ...
%     ( kappa2^(4-2*alpha) - kappa1^(4-2*alpha) );
norm_tr_a_theo = model.sigma.offset_spectrum_a_sigma * ...
    (1/(4-2*alpha)) * ...
    ( kappa2^(4-2*alpha) - kappa1^(4-2*alpha) );

% % Parseval ( * prod(model.grid.MX) )
% % and integrated over the space ( * prod(model.grid.MX) )
% norm_a_theo = prod(model.grid.MX)^2 * norm_a_theo;
%
% % Influence of discretisation
% norm_a_theo = norm_a_theo / ( prod(model.grid.MX.*model.grid.dX) /(2*pi) ) ;
%
%
% % % Influence of the complex brownian variance
% % norm_a_theo = 1/(prod(model.grid.MX))*norm_a_theo;
%
% % % Calcul of energy
% % % One has to divid by prod(model.grid.MX) because of the form of Parseval
% % % theorem for discrete Fourier transform
% % norm_a_theo = 1/prod(model.grid.MX) * norm_a_theo

end
