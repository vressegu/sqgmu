function tau = characteristic_time(model, Variance)
%% function tau = characteristic_time(model, Variance)
% Computes the characteristic time for the stochastic diffusion tensor.
%
% Written by P. DERIAN 2016-08-30.

% root mean square velocity fluctuations
% from the variance
rms1 = sqrt(mean2(Variance(:,:,1)));
rms2 = sqrt(mean2(Variance(:,:,2)));

% compute times
l1 = model.grid.dX(1);
l2 = model.grid.dX(2);
tau1 = (l1/rms1)^(2./3.)*model.advection.dt_adv^(1./3.);
tau2 = (l2/rms2)^(2./3.)*model.advection.dt_adv^(1./3.);

% return norm (mean? min?)
% tau = mean([tau1 tau2]);
tau = min([tau1 tau2]);
end