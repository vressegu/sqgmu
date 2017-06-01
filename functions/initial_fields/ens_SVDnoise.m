function buoy_ens = ens_SVDnoise(model)
%% buoy_ens = ens_SVDnoise(model)
% Generate an ensemble of low-res realizations using SVD noise.
%
% Arguments:
%   - model: the current simulation model of grid size [M, N];
% with expected model.init attributes:
%   - type_data: the name of the initial condition;
%   - P: the patch size (EVEN)
%   - N_obs: the number of pseudo-observations to generate
%   - boundary_condition: the boundary condition.
%   - scaling: the noise scaling (0=default automatic scaling).
% Output: a [M, N, 1, N_ech] array.
%
% Note: don't know why it has to be a 4-d array ... to match with velocity?
%       [TODO] ask Valentin.
%
% Written by P. DERIAN 2017-04-10

% the mean field: create a minimal model at appropriate resolution
model_tmp = fct_physical_param();
model_tmp.is_stochastic = false;
model_tmp.resolution = model.resolution;
model_tmp.init = model.init;
model_tmp.init.type_rand = 'None';
model_tmp.advection.N_ech = 1;
[fft_buoy_mean, ~] = fct_buoyancy_init_ens(model_tmp);
buoy_mean = real(ifft2(fft_buoy_mean));

% create noise instance
noise_generator = SVDnoise3D(model.grid.MX, model.init.P, ...
                             model.init.boundary_condition);
% compute the noise basis from this mean field                         
scaling = model.init.scaling;
if scaling==0
    scaling = noise_generator.default_amplitude_scaling();
end
[Phi, Sigma] = noise_generator.scalar_noise_basis(buoy_mean, model.init.N_obs, ...
                                                  scaling);
gamma = randn(numel(Sigma), model.advection.N_ech);
buoy_ens = Phi*bsxfun(@times, gamma, Sigma);

% restore mean
buoy_ens = bsxfun(@plus, buoy_ens, buoy_mean(:));

% reshape
buoy_ens = reshape(buoy_ens, [size(buoy_mean), 1, model.advection.N_ech]);
