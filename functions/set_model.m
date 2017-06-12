function model = set_model()
%% Set the model parameters
%
% NB: [WIP] tags stand for Work In Progress, beware when using these
% options.
%
% NB: reproductability is NOT GUARANTEED at the moment for stochastic
% simulations due to the different random generators on each worker.
%
% Added by P. DERIAN 2016-10-12.
% Modified by P. DERIAN many times :)
% Modified by P. DERIAN 2017-04-10:
%   - added 'Scramble" init and the "init" sub-structure in "model".
% Modified by P. DERIAN 2017-05-03: 
%   - added the "type_rand" parameter;
%   - split the initial condition form its optional randomization.
% TODO: 
%   - sometimes I feel like N_ech should not be a sub-attribute of
%       model.advection, but rather of model.init or at model root...
%   - seems like the parfor in fct_fft_advection_sto() messes up with the
%   random seed. fix it!

%% Main parameters

% Deterministic or random model
%------------------------------
stochastic_simulation = 1;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

% Duration of the simulation (in seconds)
%----------------------------------------
advection_duration = 3600*24*30; % in [s] (last number is days)

% Number of realizations in the ensemble
%---------------------------------------
N_ech = 12;
% NB: 
% - with deterministic simulation (stochastic_simulation = 1) and
%   ensemble run (N_ech > 1), make sure set a randomization method
%   (type_rand below), otherwise all particles are the same.
% - N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128.

% Type of initial conditions
%---------------------------
type_data = 'Vortices2';
% - 'Vortices' : 2 large anticyclones and 2 large cyclones
%    (used in "Geophysical flow under location uncertainty", Resseguier V.,
%     Memin E., Chapron B.)
% - 'Vortices2' : same as 'Vortices' but properly periodized.
% - 'Perturbed_vortices' : Same flow with slight small-scale modifications
%    (used in "Chaotic transitions and location uncertainty in geophysical
%    fluid flows", Resseguier V., Memin E., Chapron B.)
% - 'Spectrum' : Gaussian random field with a spectrum slope defined by
%    the variable slop_b_ini (default value  = -5/3)
% - [WIP] 'PertubVortices2': NOT COMPATIBLE with type_rand={'Scramble', 'ScrambleSc'}

% Type of initialization randomization
%-------------------------------------------
type_rand = 'None';
% - 'None' or '': no randomization.
% - 'SVDnoise': use the SVD-based, same as for the 'SVDfull' stochastic model.
% - 'Scramble': randomize from a high-resolution observation of the initial
%    condition (true variance).
% - 'ScrambleSc': same as 'Scramble' but allows to scale the noise variance.
% - [WIP] 'Spectrum': add a small-scale spectral noise of prescribed slope.

% Type of stochastic noise (stochastic_simulation=1)
%---------------------------------------------------
type_noise = 'SVDfull';
% - 'None' or '': no noise, purely deterministic.
% - 'Spectrum': isotropic, homogeneous small-scale spectral noise
% - 'SVDfull': noise basis learned from pseudo-observations, built by SVD.

% [WIP] Type of forcing
%----------------------
type_forcing = 'None';
% NB: this has not been tested, DO NOT USE unless you know what you're
%    doing.
% - 'None': no forcing on the large-scale velocity w;
% - 'Jet': a constant jet.

% Resolution and grid
%--------------------
resolution = 128;
% The number of grid point is resolution^2
% It has to be an even integer
dealias_method = 'exp';
% [WIP] Method for mandatory de-aliasing of non-linear terms in
% pseudospectral codes (advection and non-homogeneous stochastic diffusion) 
% - 'lowpass': same as in SQGMU 1;
% - '2/3': the classical 2/3 rule (NB: produces Gibb's oscillations);
% - 'exp': high-order exponential filter (Constantin et al., J. Sci.
%   Comput. (2012)).

% Random generator seed
%----------------------
% if <1, ignored.
% [TODO] deal with parfor & workers random generators
seed = 19860302; % positive integer.

%% Optional parameters

% Output controls
%----------------
verbose = 1; % verbose level
plot_results = 1; % toggle display
folder_output = '/Users/pderian/Documents/Data/simu_SQG/TMP'; % root directory for output data
plot_moments = 1; % plot one-point one-time moments each day (requires plot_results)
image_format = 'png'; % format of output images - see "print formattype".
                      % Use 'epsc' for color EPS files.
image_dpi = 150; % resolution (dot per inch) of output images.

% Parameters for the chosen stochastic noise model
%--------------------------------------------------
sigma = struct();
switch type_noise
    case 'Spectrum'
        % Isotropic homogeneous spectral noise.
        %-----------------------------------------
        k_c = 1./(3e2); %1/(300 meters)
        % a_H is calculated latter in the code using
        % a_H = 2 * f_0 / k_c^2
        % where f_0 is the Coriolis frequency
        % Spectrum slope of sigma dBt
        slope_sigma = -5./3.;
        % Rate between the smallest and the largest wave number of sigma dBt
        kappamin_on_kappamax = 1./2.;
        %-----------------------------------------
        % Gather sigma parameters in one structure
        % NB: do NOT edit values in the remaining of this section
        sigma.k_c = k_c;
        sigma.slope_sigma = slope_sigma;
        sigma.kappamin_on_kappamax = kappamin_on_kappamax;
        sigma.long_name = 'Spectral noise';
    case 'SVDfull'
        % Model using a noise basis built in-line from pseudo-observations by SVD.
        %-----------------------------------------
        k_c = 1./(3e2); %1/(300 meters) [TODO] remove?
        P = 3; %observation patch side size - ODD value
        N_obs = 20; %number of pseudo-observations to generate
        boundary_condition = 'circular'; %patch boundary condition
        % project the complete velocity field to divergence-free when
        % computing the advection
        divfree_projection = false; %[WIP] not fully tested
        %-----------------------------------------
        % Gather sigma parameters in one structure
        % NB: do NOT edit values in the remaining of this section
        sigma.k_c = k_c; % [TODO] remove?
        sigma.P = P;
        sigma.N_obs = N_obs;
        sigma.boundary_condition = boundary_condition;
        sigma.divfree_projection = divfree_projection;
        sigma.long_name = 'Pseudo-observations noise';
    case {'None', ''}
        % No noise, deterministic simulation
        %-----------------------------------------
        % If the simulation is deterministic, a_H = 0 and only one simulation
        % is performed
        k_c = inf; % And then a_H = 0
        % [TODO] is k_c necessary with deterministic? or was it left here just
        % to test Stoch or Det? ask Val and maybe remove it.
        %-----------------------------------------
        % Gather sigma parameters in one structure
        % Note: do NOT edit values in the remaining of this section
        sigma.long_name = 'None';
        sigma.k_c = k_c;
    % Model not defined, abort.
    otherwise
        error('SQGMU:main:ValueError', '"%s" is not a valid noise model', type_noise);
end
% also remember the noise type
sigma.type_noise = type_noise;

% Parameters for the initial condition and its randomization
%-----------------------------------------------------------
init = struct();
% the data
switch type_data
    case {'Vortices', 'Vortices2', 'Pertubed_vortices'}
        % Nothing to do here
    case 'PerturbVortices2'
        init.delta_err = 0.1; % the std, in [grid point], of vortex center error on each direction
        init.N_avg = 20; % the number of samples to be averaged
    case 'Spectrum'
        init.slope_b_ini = -5./3.; % spectrum slope of the initial condition
    otherwise
        error('sqgmu:set_model:InvalidParameter', 'the initial condition "%s" is unknown', type_data);
end
init.type_data = type_data;
% the randomization
switch type_rand
    case {'None', ''}
        % Nothing to do here
    case 'Scramble'
        init.resolution_ratio = 4; %ratio between high-res and simulation resolution
        init.circshift = false; %shift the observation patch across the periodic boundary
    case 'ScrambleSc'
        init.resolution_ratio = 4; %ratio between high-res and simulation resolution
        init.circshift = false; %shift the observation patch across the periodic boundary
        init.scaling = 0.50; %scalign factor for the variance
    case 'SVDnoise'
        init.P = 3; %patch size (ODD)
        init.N_obs = 50; %number of pseudo-observations to generate
        init.boundary_condition = 'circular'; %patch boundary condition
        init.scaling = 0; %if <=0, use the SVDNoise default scaling.
    case 'Spectrum'
        init.slope_b_ini = -5./3.; %[TODO] could be conflicting with that of type_data
        init.k_min = 2e-4; %[TODO] as a ratio kmin/kmax
        init.scaling = 1e-4; %[TODO] make it simpler
end;
init.type_rand = type_rand;

% Parameters for the forcing
%---------------------------
forcing = struct();
switch type_forcing
    case 'None'
        forcing.amplitude = 0.;
        % do nothing
    case 'Jet'
        forcing.std = 2e5; % std of the gaussian curve for jet width
        forcing.amplitude = .5; % amplitude (and sign) of the jet
end
forcing.type_forcing = type_forcing;

%% Gather all parameters

% Physical parameters
model = fct_physical_param();

% Other parameters
model.is_stochastic = stochastic_simulation; % whether stochastic or deterministic simulation
model.sigma = sigma; % stochastic parameters
model.init = init; % the initial condition parameters
model.type_data = model.init.type_data; % left for backward-compatibility [TODO] remove?
model.resolution = resolution; %grid resolution
model.grid.dealias_method = dealias_method; %de-aliasing method
model.advection.N_ech = N_ech; %number of simultaneous runs
model.advection.advection_duration = advection_duration; %simulation time
model.advection.forcing = forcing; %the forcing of the velocity
model.verbose = verbose; %verbose level
model.output.folder_root = folder_output; %root folder for output
model.output.plot_results = plot_results; %toggle plotting
model.output.plot_moments = plot_moments; %plotting options
model.output.image_format = image_format; %output graphics format
model.output.image_dpi = image_dpi; %output graphics resolution

%% Random generator
% if no seed was specified, generate from current time.
if seed<1
    seed = rem(ceil(now*1e6), 2^32); %[TODO] likely not so good?
end
rng(seed);
model.seed = seed;
%[TODO] set seed to the different workers...

return
