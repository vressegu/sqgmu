function model = set_model()
%% Set the model parameters
%
% Added by P. DERIAN 2016-10-12.

%% Main parameters to choose

% Deterministic or random model
stochastic_simulation = 1;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

% Duration of the simulation (in seconds)
advection_duration = 3600*24*80; % last number is days

% Number of realizations in the ensemble
N_ech = 200;
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations with type_noise='None' )

% Type of initial condtions 
type_data = 'Vortices';
% 'Vortices' : 2 large anticyclones and 2 large cyclones
%   (used in "Geophysical flow under location uncertainty", Resseguier V.,
%    Memin E., Chapron B.)
% 'Perturbed_vortices' : Same flow with slight small-scale modifications
%   (used in "Chaotic transitions and location uncertainty in geophysical 
%    fluid flows", Resseguier V., Memin E., Chapron B.)
% 'Spectrum' : Gaussian random field with a spectrum slope defined by 
%   the variable slop_b_ini (default value  = -5/3)

% Type of random noise:
%    - for the stochastic model (stochastic_simulation=1)
%    - or to randomize initial condition (stochastic_simulation=0, N_ech>1)
type_noise = 'SVDfull';
% 'Spectrum' : isotropic, homogeneous small-scale spectral noise
% 'SVD' : noise basis learned from pseudo-observations, built by SVD.
% 'None' : no noise, purely deterministic. N_ech set to 1.

% [WIP] Type of forcing
type_forcing = 'None';
% 'None': no forcing on the large-scale velocity w;
% 'Jet': a constant jet.

% Resolution
resolution = 128;
% The number of grid point is resolution^2
% It has to be an even integer

% Random generator seed
% if <1, ignored.
% [TODO] check si pas fixe ailleurs
seed = 19860302; % positive integer.

%% Optional parameters

% Output controls
verbose = 1; % verbose level
plot_results = 1; % toggle display
folder_output = '/Users/pderian/Documents/Data/simu_SQG/test_parr'; % root directory for output data
plot_moments = 1; % plot one-point one-time moments each day (requires plot_results)
image_format = 'png'; % format of output images - see "print formattype".
                      % Use 'epsc' for color EPS files.
image_dpi = 150; % resolution (dot per inch) of output images.

% Switch according to the chosen noise model
sigma = struct();
switch type_noise
    % Isotropic homogeneous spectral noise.
    case 'Spectrum'
        k_c = 1./(3e2); % 1/(300 meters)
        % a_H is calculated latter in the code using 
        % a_H = 2 * f_0 / k_c^2
        % where f_0 is the Coriolis frequency
        % Spectrum slope of sigma dBt
        slope_sigma = -5./3.;
        % Rate between the smallest and the largest wave number of sigma dBt
        kappamin_on_kappamax = 1./2.;
        % Gather sigma parameters in one structure
        % Note: do NOT edit values in the remaining of this section
        sigma.k_c = k_c;
        sigma.slope_sigma = slope_sigma;
        sigma.kappamin_on_kappamax = kappamin_on_kappamax;
        sigma.long_name = 'Spectral noise';
    % Model using a noise basis built in-line from pseudo-observations by SVD.
    case {'SVD', 'SVDfull'}
        k_c = 1./(3e2); % 1/(300 meters) [TODO] remove?
        % patch side size
        P = 5; % EVEN for SVD, ODD for SVDfull
        % number of pseudo-observations to generate
        N_obs = 20;      
        % interpolation method for variance data
        interp_method = 'bicubic'; % for 'SVD'
        boundary_condition = 'circular'; % for 'SVDfull'
        % project the complete velocity field to divergence-free when
        % computing the advection
        divfree_projection = false;
        % Gather sigma parameters in one structure
        % Note: do NOT edit values in the remaining of this section
        sigma.k_c = k_c; % [TODO] remove?
        sigma.P = P;
        sigma.N_obs = N_obs;
        sigma.interp_method = interp_method;
        sigma.boundary_condition = boundary_condition;
        sigma.divfree_projection = divfree_projection;
        sigma.long_name = 'Pseudo-observations noise';
    % Noise previously learnt on fine-resolution data (also with SVD) 
    % [TODO] fix
    case 'Learning_SVD'
         k_c = 1./(3e2); % 1/(300 meters) [TODO] remove?
        % data location
        observation_dir = '/Volumes/Calcul/pderian/SQGMU/usual_SQG/Vortices_256';
        observation_pattern = 'var*.mat';
        % project the complete velocity field to divergence-free when
        % computing the advection
        divfree_projection = false;
        % Gather sigma parameters in one structure
        % Note: do NOT edit values in the remaining of this section
        sigma.k_c = k_c; % [TODO] remove?
        sigma.observation_dir = observation_dir;
        sigma.observation_pattern = observation_pattern;
        sigma.divfree_projection = divfree_projection;
        sigma.long_name = 'Noise learnt on fine-resolution data with SVD';
    % No noise
    case 'None'
        % If the simulation is deterministic, a_H = 0 and only one simulation
        % is performed
        k_c = inf; % And then a_H = 0
        % [TODO] is k_c necessary with deterministic? or was it left here just
        % to test Stoch or Det? ask Val and maybe remove it.
        % Gather sigma parameters in one structure
        % Note: do NOT edit values in the remaining of this section
        sigma.type_noise = 'None';
        sigma.long_name = 'None';
        sigma.k_c = k_c;
        % Force some other parameters
        N_ech = 1;
        plot_moments = false; 
    % Model not defined, abort.
    otherwise
        error('SQGMU:main:ValueError', '"%s" is not a valid noise model', type_noise);
end
% also remember the noise type
sigma.type_noise = type_noise;    

% Parameters for the forcing
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

% Spectrum slope of the initial condition (if type_data = 'Spectrum' )
slope_b_ini = -5./3.; 

%% Gather all parameters

% Physical parameters
model = fct_physical_param();

% Other parameters
model.is_stochastic = stochastic_simulation; % whether stochastic or deterministic simulation
model.sigma = sigma; % stochastic parameters
model.type_data = type_data; % initial condition
if strcmp(type_data,'Spectrum') % initial condition
    model.slope_b_ini = slope_b_ini;
end
model.resolution = resolution; %grid resolution
model.advection.N_ech = N_ech; % number of simultaneous runs
model.advection.advection_duration = advection_duration; % simulation time
model.advection.forcing = forcing; % the forcing of the velocity
model.verbose = verbose; % verbose level
model.output.folder_root = folder_output; % root folder for output
model.output.plot_results = plot_results; % toggle plotting 
model.output.plot_moments = plot_moments; % plotting options
model.output.image_format = image_format; % output graphics format
model.output.image_dpi = image_dpi; %output graphics resolution

%% Random generator
% if no seed was specified, generate from current time.
if seed<1
    seed = ceil(now*1e6);
end
rng(seed);
model.seed = seed;

return 