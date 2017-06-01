function buoy_ens = ens_Spectrum(model)
%% buoy_ens = ens_Spectrum(model)
% Generate an ensemble of low-res realizations by adding high-freq  spectral noise.
%
% Arguments:
%   - model: the current simulation model of grid size [M, N];
% with expected model.init attributes:
%   - type_data: the name of the initial condition;
%   - slope_b_ini: the spectrum slope
%   - k_min: the min wavenumber
%   - scaling: a scaling factor
%
% Output: a [M, N, 1, N_ech] array.
%
% Note: don't know why it has to be a 4-d array ... to match with velocity?
%       [TODO] ask Valentin.
%
% Written by P. DERIAN 2017-05-04

% the mean field: create a minimal model at appropriate resolution
model_tmp = fct_physical_param();
model_tmp.is_stochastic = false;
model_tmp.resolution = model.resolution;
model_tmp.init.type_data = model.init.type_data;
model_tmp.init.type_rand = 'None';
model_tmp.advection.N_ech = 1;
[fft_buoy_mean, ~] = fct_buoyancy_init_ens(model_tmp);

% replicate
buoy_ens = repmat(real(ifft2(fft_buoy_mean)), [1, 1, 1, model.advection.N_ech]);

% add spectrum
for n=1:model.advection.N_ech
    buoy_ens(:,:,:,n) = buoy_ens(:,:,:,n) + init_Spectrum(model, model.init.k_min, model.init.scaling); % with gap
end

