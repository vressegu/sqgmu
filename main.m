%%%%%%%%%%%%%%%%%%%%
%%% Main 
%%%%%%%%%%%%%%%%%%%%
init;

%% Main parameters to choose

% Type of dynamics
dynamics = 'SQG';
% dynamics = '2D';

% Deterministic or random model
stochastic_simulation = false;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

% Duration of the simulation (in seconds)
advection_duration = 3600*24*30;
% advection_duration = 3600*24*1000;
% % advection_duration = 3600*24*20; % 20 days

% Number of realizations in the ensemble
N_ech=1;
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )

% Type of initial condtions 
type_data ='Vortices';
% 'Vortices' : 2 large anticyclones and 2 large cyclones
%   (used in "Geophysical flow under location uncertainty", Resseguier V.,
%    Memin E., Chapron B.)
% 'Vortices2' : same as 'Vortices' but properly periodized (by Pierre Derian).
% 'Perturbed_vortices' : Same flow with slight small-scale modifications
%   (used in "Chaotic transitions and location uncertainty in geophysical 
%    fluid flows", Resseguier V., Memin E., Chapron B.)
% 'Spectrum' : Gaussian random field with a spectrum slope deined by 
%   the variable slop_b_ini (default value  = -5/3)
% 'Zero' : Field equal to zero everywhere

% Resolution
resolution = 128;
%resolution = 512;
%resolution = 128;
%resolution = 1024;

% The number of grid point is resolution^2
% It has to be an even integer

% Boundaries conditions
dirichlet = false;

% Forcing

% Forcing or not
forcing = false;
% If yes, there is a forcing
% F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y)
% % If yes, there is an additionnal velocity V = (0 Vy) 
% % with Vy = ampli_forcing * odg_b *  sin( 2 freq_f pi y/L_y)

% Type de forcing
% forcing_type = 'Kolmogorov';
forcing_type = 'Spring';

% Amplitude of the forcing
ampli_forcing = 10;
% ampli_forcing = 1;

% Frequency of the forcing
freq_f = [3 2];
% freq_f = [0 1];

% Viscosity
Lap_visco.bool = false;
% HV.bool = false;

% % Smagorinsky-like viscosity
% Smag.bool = false;
% % HV.bool = false;

% Hyper-viscosity
HV.bool = true;
% HV.bool = false;

% Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
Smag.bool = true;

% For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
% Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) ) 
% and the targeted diffusion scale
Smag.kappamax_on_kappad = 0.8;
% Smag.kappad_on_kappamax = 1/2;
warning('This value needed to be tuned?')

%% Optional parameters

% Choose to plot one-point one-time moments each day
plot_moments = true;

% Plot dissipations terms
plot_dissip = true;

% Begin simulation from a precomputed field?
use_save = false;
% In this case, which day should be used as initialisation
day_save = 96;

% Variance tensor a_H
if stochastic_simulation
    switch dynamics
        case 'SQG'
            k_c = 1/(3e2); % 1/(300 meters)
        case '2D'
            error(...
      'The turbulence 2D is not stable under the action of noise');
%             k_c = 1/(eps); % 1/(100 meters)
        otherwise
            error('Unknown type of dynamics');
    end
    % a_H is calculated latter in the code using
    % a_H = 2 * f_0 / k_c^2
    % where f_0 is the Corilis frequency
else
    % If the simulation is deterministic, a_H = 0 and only one simulation
    % is performed
    k_c = inf; % And then a_H = 0
    N_ech=1;
    plot_moments = false; 
end

% Spectrum slope of sigma dBt
switch dynamics
    case 'SQG'
        slope_sigma = - 5/3; 
    case '2D'
        slope_sigma = - 3;
    otherwise
        error('Unknown type of dynamics');
end

% Rate between the smallest and the largest wave number of sigma dBt
kappamin_on_kappamax = 1/2;

% Spectrum slope of the initial condition (if type_data = 'Spectrum' )
switch dynamics
    case 'SQG'
        slope_b_ini = - 5/3;
    case '2D'
        slope_b_ini = - 3;
    otherwise
        error('Unknown type of dynamics');
end

% Physical parameters
model = fct_physical_param(dynamics);

% Gather parameters in the structure model
model.sigma.slope_sigma = slope_sigma;
model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
if strcmp(type_data,'Spectrum')
    model.slope_b_ini = slope_b_ini;
end
model.dynamics=dynamics;
model.type_data=type_data;
model.advection.N_ech=N_ech;
model.sigma.k_c = k_c;
model.advection.advection_duration=advection_duration;
model.advection.plot_dissip = plot_dissip;
model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
model.advection.forcing.ampli_forcing = ampli_forcing;
model.advection.forcing.freq_f = freq_f;
model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag;
model.advection.use_save = use_save;
model.advection.day_save = day_save;

%% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto(model, fft_buoy);
