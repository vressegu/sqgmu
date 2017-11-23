function main(stochastic_simulation,type_data,resolution,forcing, ...
    Lap_visco,HV,Smag)
% Main function to Launch thte code
%

%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    init;
end

%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
%dynamics = '2D';

% Type of spectrum for sigma dBt
%type_spectrum = 'Band_Pass_w_Slope'; % as in GAFD part II
%type_spectrum = 'Low_Pass_w_Slope';
% Spectrum cst for k<km ans slope for k>km
% type_spectrum = 'Low_Pass_streamFct_w_Slope';
% Matern covariance for the streamfunction
% spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
% ~ k2 for k<km ans slope for k>km
% type_spectrum = 'BB';
% type_spectrum = 'Bidouille';
type_spectrum = 'SelfSim_from_LS';
%  Sigma computed from self similarities from the large scales
sigma.type_spectrum = type_spectrum;


if nargin == 0
    % Deterministic or random model
    stochastic_simulation = true;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
    
    % Smagorinsky-like ocntrol of dissipation
    sigma.Smag.bool = false;
    
    %     % Sigma computed from self similarities from the large scales
    %     sigma.SelfSim_from_LS.bool = true;
    
    %     if sigma.SelfSim_from_LS.bool
    %         % Sigma computed from a energy of absolute diffusivity spectrum
    %         % sigma.SelfSim_from_LS.spectrum = 'energy';
    %         sigma.SelfSim_from_LS.spectrum = 'abs_diff';
    %     end
    
    if strcmp(type_spectrum,'SelfSim_from_LS')
        % Heterrogeenosu energy flux epsilon
        sigma.hetero_energy_flux = false;
        
        % Modulation by local V L (estimated from the velocity and from 
        % thegradient of the velocity)
        sigma.hetero_modulation = false;
        
        % Modulation by local V^2 
        sigma.hetero_modulation_V2 = true;
        
        %     %if strcmp(type_spectrum,'SelfSim_from_LS')
        %     if sigma.hetero_modulation & strcmp(type_spectrum,'SelfSim_from_LS')
        if sigma.hetero_modulation | sigma.hetero_energy_flux ...
                | sigma.hetero_modulation_V2
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            % Smag.dealias_ratio_mask_LS = 1/16;
            Smag.dealias_ratio_mask_LS = 1/8;
            
        end
    end
    
    % Force sigma to be diveregence free
    sigma.proj_free_div = true;
    
    if (sigma.Smag.bool + sigma.hetero_modulation + ...
            sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ) > 1
        error('These two parametrizations cannot be combined');
    end        
    
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if sigma.Smag.bool
        % Smagorinsky energy budget (dissipation epsilon)
        % without taking into account the noise intake
        sigma.Smag.epsi_without_noise = false;
        
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1;
        % Smag.dealias_ratio_mask_LS = 1/8;
% %         Smag.dealias_ratio_mask_LS = 1/64;
% % %         Smag.dealias_ratio_mask_LS = 1/8;
% % % %         Smag.dealias_ratio_mask_LS = 1/128;
% % % %         %         sigma.Smag.dealias_ratio_mask_LS = 1/128;
% % % %         % %         sigma.Smag.dealias_ratio_mask_LS = 1/8;
        
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
% %        sigma.Smag.kappamax_on_kappad = 2;
%        sigma.Smag.kappamax_on_kappad = 1;
        sigma.Smag.kappamax_on_kappad = 0.5;
        
%         % Factor in front of the additional constant dissipation
%         % Set to 0 for no additional constant dissipation
%         sigma.Smag.weight_cst_dissip = 0;
        
        % Rate between the smallest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMinUnresolved_on_kappaShanon = 1;
        
        % Rate between the largest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMaxUnresolved_on_kappaShanon = 8;
        
        % Heterogeneity of the noise
        sigma.Smag.SS_vel_homo = true;
        
    end
end

% Number of realizations in the ensemble
N_ech=1;
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )

% Duration of the simulation (in seconds)
advection_duration = 3600*24*30;
% advection_duration = 3600*24*1000;
% % advection_duration = 3600*24*20; % 20 days

if nargin == 0
    % Type of initial condtions
    type_data = 'Constantin_case2';
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
    % 'Constantin_case1'
    % 'Constantin_case2'
    
    % Resolution
    %resolution = 64;
    resolution = 128;
    %resolution = 256;
    %resolution = 512;
    %resolution = 1024;
    % resolution = 2048;
    
    % The number of grid point is resolution^2
    % It has to be an even integer
    
    % Forcing
    
    % Forcing or not
    forcing = false;
    % If yes, there is a forcing
    % F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y)
    % % If yes, there is an additionnal velocity V = (0 Vy)
    % % with Vy = ampli_forcing * odg_b *  sin( 2 freq_f pi y/L_y)
end

% Type de forcing
% forcing_type = 'Kolmogorov';
forcing_type = 'Spring';

% Amplitude of the forcing
ampli_forcing = 10;
% ampli_forcing = 1;

% Frequency of the forcing
freq_f = [3 2];
% freq_f = [0 1];

if nargin == 0
    % Viscosity
    Lap_visco.bool = false;
    
    % % Smagorinsky-like viscosity
    % Smag.bool = false;
    % % HV.bool = false;
    
    % Hyper-viscosity
    HV.bool = true;
    
    % Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
    Smag.bool = false;
        
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if Smag.bool
        if Lap_visco.bool
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            %     dealias_ratio_mask_LS = 1/8;
            %     %dealias_ratio_mask_LS = 1/2;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            %     %     Smag.kappamax_on_kappad = 2;
            %     % % Smag.kappamax_on_kappad = 1.5; % Better, still small oscillations or just pixels?
            %     % %    %  Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
            %     Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
            %     %  d'aliasing
            Smag.kappamax_on_kappad = 1; % Stable mais petit artefact
            %  d'aliasing
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            %    HV.weight_cst_dissip = 1/10;
            %     HV.weight_cst_dissip = 1/10; % no aliasing
            Smag.weight_cst_dissip = 0;
        elseif HV.bool
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            Smag.kappamax_on_kappad = 1.1;% still small oscillations or just pixels?
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
            %     % HV.weight_cst_dissip = 1/3; % still a bit of (stable) aliasing
            %     HV.weight_cst_dissip = 1/3;
            %     % HV.weight_cst_dissip = 0;
            
            Smag.weight_cst_dissip = 1/1;
            % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
        end
    else
        Smag.kappamax_on_kappad = 0;
    end
    %     % Smag.kappamax_on_kappad = 1.1;
    %     % % Smag.kappad_on_kappamax = 1/2;
    %     if Smag.bool
    %         warning('This value needed to be tuned?')
    %     end
end

%% Optional parameters
% Choose to plot
if nargin > 0
    plots_bool = false;
else
    plots_bool = true;
end

% Compute velocity covariance and absolute diffusivity
cov_and_abs_diff = false;

% Choose to plot one-point one-time moments each day
plot_moments = true;

% Choose to plot the dissipation by scale
plot_epsilon_k = false;
if sigma.hetero_energy_flux
    plot_epsilon_k = true;
end

% Plot dissipations terms
plot_dissip = true;

% Begin simulation from a precomputed field?
use_save = false;
% In this case, which day should be used as initialisation
day_save = 7;

dealias_method = 'exp';
% [WIP] Method for mandatory de-aliasing of non-linear terms in
% pseudospectral codes (advection and non-homogeneous stochastic diffusion)
% - 'lowpass': same as in SQGMU 1;
% - '2/3': the classical 2/3 rule (NB: produces Gibb's oscillations);
% - 'exp': high-order exponential filter (Constantin et al., J. Sci.
%   Comput. (2012)).

% Boundaries conditions
dirichlet = false;

% Variance tensor a_H
if stochastic_simulation
    switch dynamics
        case 'SQG'
            sigma.k_c = 1/(3e2); % 1/(300 meters)
        case '2D'
            if strcmp(type_spectrum , 'SelfSim_from_LS')
                sigma.k_c = 0;
            else
                error(...
                    'The turbulence 2D is not stable under the action of noise');
                %             k_c = 1/(eps); % 1/(100 meters)
            end
        otherwise
            error('Unknown type of dynamics');
    end
    % a_H is calculated latter in the code using
    % a_H = 2 * f_0 / k_c^2
    % where f_0 is the Corilis frequency
else
    % If the simulation is deterministic, a_H = 0 and only one simulation
    % is performed
    sigma.k_c = inf; % And then a_H = 0
    N_ech=1;
    plot_moments = false;
end

% Spectrum slope of sigma dBt
switch dynamics
    case 'SQG'
        sigma.slope_sigma = - 5/3;
    case '2D'
        sigma.slope_sigma = - 3;
    otherwise
        error('Unknown type of dynamics');
end
if  strcmp(sigma.type_spectrum,'BB')
    sigma.slope_sigma = 0;
end

% Rate between the smallest and the largest wave number of sigma dBt
if strcmp(type_spectrum , 'SelfSim_from_LS')
    sigma.kappamin_on_kappamax = 1/2;
    % sigma.kappamin_on_kappamax = 1/8;
    sigma.kappaLS_on_kappamax = 1/8;
else
    %kappamin_on_kappamax = 1/32;
    sigma.kappamin_on_kappamax = 1/2;
    % sigma.kappamin_on_kappamax = 1/128;
    %         sigma.slope_sigma = - 5;
    % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
end

% Rate between the largest wave number of sigma dBt and the largest wave
% number of the simulation
sigma.kappamax_on_kappaShanon = 1;

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
model.sigma = sigma;
eval(['model.sigma.fct_tr_a = @(m,k1,k2,alpha) fct_norm_tr_a_theo_' ...
    model.sigma.type_spectrum '(m,k1,k2,alpha);']);
% model.sigma.slope_sigma = slope_sigma;
% model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
if strcmp(type_data,'Spectrum')
    model.slope_b_ini = slope_b_ini;
end
model.dynamics=dynamics;
model.type_data=type_data;
model.advection.N_ech=N_ech;
% model.sigma.k_c = k_c;
model.advection.advection_duration=advection_duration;
model.advection.plot_epsilon_k = plot_epsilon_k;
model.advection.plot_dissip = plot_dissip;
model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
model.advection.forcing.ampli_forcing = ampli_forcing;
model.advection.forcing.freq_f = freq_f;
model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.cov_and_abs_diff = cov_and_abs_diff;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag;
model.advection.use_save = use_save;
model.advection.day_save = day_save;
model.grid.dealias_method = dealias_method; %de-aliasing method
%model.Smag.dealias_ratio_mask_LS = dealias_ratio_mask_LS;
model.plots = plots_bool;

%% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto(model, fft_buoy);
