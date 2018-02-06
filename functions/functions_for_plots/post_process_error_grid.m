function error_vs_t = post_process_error_grid(stochastic_simulation,...
    type_data,resolution,resolution_HR,forcing,sigma,Lap_visco,HV,Smag,...
    N_ech,first_day)
%Lap_visco,HV,Smag,day_choose)
% plot the same thing that fct_fft_advection_sto
%


if nargin == 0
    init;
end
%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
%dynamics = '2D';

plot_random_IC = false;
random_IC_large = false


if nargin == 0
    
    % Deterministic or random model
    stochastic_simulation = true;
    sigma.sto = stochastic_simulation;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
    
    if sigma.sto
        % Type of spectrum for sigma dBt
        % sigma.type_spectrum = 'Band_Pass_w_Slope' % as in GAFD part II
        % sigma.type_spectrum = 'Low_Pass_w_Slope';
        % Spectrum cst for k<km ans slope for k>km
        % sigma.type_spectrum = 'Low_Pass_streamFct_w_Slope';
        % Matern covariance for the streamfunction
        % spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
        % ~ k2 for k<km ans slope for k>km
        % type_spectrum = 'BB';
        % type_spectrum = 'Bidouille';
        % sigma.type_spectrum = 'SelfSim_from_LS'
        sigma.type_spectrum = 'EOF'
        %  Sigma computed from self similarities from the large scales
        % sigma.type_spectrum = type_spectrum;
        
        % Homogeneous dissipation associated with the spectrum slope
        sigma.assoc_diff = false;
        
        % Smagorinsky-like control of dissipation
        sigma.Smag.bool = false;
        
        %     % Sigma computed from self similarities from the large scales
        %     sigma.SelfSim_from_LS.bool = true;
        
        %     if sigma.SelfSim_from_LS.bool
        %         % Sigma computed from a energy of absolute diffusivity spectrum
        %         % sigma.SelfSim_from_LS.spectrum = 'energy';
        %         sigma.SelfSim_from_LS.spectrum = 'abs_diff';
        %     end
        
        % if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        % Heterrogeenosu energy flux epsilon
        sigma.hetero_energy_flux = false;
        
        % Modulation by local V L (estimated from the velocity and from
        % thegradient of the velocity)
        sigma.hetero_modulation = false;
        
        % Modulation by local V^2
        sigma.hetero_modulation_V2 = false;
        
        %     %if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        %     if sigma.hetero_modulation & strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        if sigma.hetero_modulation | sigma.hetero_energy_flux ...
                | sigma.hetero_modulation_V2 || strcmp(sigma.type_spectrum,'EOF')
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            % Smag.dealias_ratio_mask_LS = 1/4;
            
        end
        % end
        
        % Force sigma to be diveregence free
        sigma.proj_free_div = true;
        
        if ( (sigma.Smag.bool + sigma.hetero_modulation + ...
                sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ) > 1 ) ...
                || ( (sigma.Smag.bool + sigma.assoc_diff ) > 1 )
            error('These parametrizations cannot be combined');
        end
        
        if sigma.Smag.bool || sigma.assoc_diff
            % Rate between the smallest wave number of the spatially-unresolved
            % (not simulated) component of sigma dBt and the largest wave
            % number of the simulation
            sigma.kappaMinUnresolved_on_kappaShanon = 1;
            
            % Rate between the largest wave number of the spatially-unresolved
            % (not simulated) component of sigma dBt and the largest wave
            % number of the simulation
            sigma.kappaMaxUnresolved_on_kappaShanon = 8;
            
        end
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        if sigma.Smag.bool
            % Smagorinsky energy budget (dissipation epsilon)
            % without taking into account the noise intake
            sigma.Smag.epsi_without_noise = false;
            
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            % Smag.dealias_ratio_mask_LS = 1;
            % Smag.dealias_ratio_mask_LS = 1/8;
            % Smag.dealias_ratio_mask_LS = 1/4;
            %Smag.dealias_ratio_mask_LS = 1/2;
            Smag.dealias_ratio_mask_LS = 1;
            
            %         % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            %         % and the targeted diffusion scale
            %         % %        sigma.Smag.kappamax_on_kappad = 2;
            %         % sigma.Smag.kappamax_on_kappad = 1;
            
            % sigma.Smag.kappamax_on_kappad = 0.5; % (better(?))
            % sigma.Smag.kappamax_on_kappad = 1 / 4;
            sigma.Smag.kappamax_on_kappad = 1 / ...
                sigma.kappaMaxUnresolved_on_kappaShanon;
            
            %         % Factor in front of the additional constant dissipation
            %         % Set to 0 for no additional constant dissipation
            %         sigma.Smag.weight_cst_dissip = 0;
            
            % Heterogeneity of the noise
            sigma.Smag.SS_vel_homo = false;
            
        end
        
        % Desactivate the noise
        sigma.no_noise = false;
        if sigma.no_noise
            warning('There is no noise here');
        end
    end
end

% Number of realizations in the ensemble
if nargin == 0
    % N_ech=1;
    N_ech=200;
    % else
    %     N_ech=1;
end
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )


% Duration of the simulation (in seconds)
advection_duration = 3600*24*30;
%advection_duration = 3600*24*1000;
% % advection_duration = 3600*24*20; % 20 days

if nargin == 0
    
    first_day = 0
    % first_day = 101
    
    % Type of initial condtions
    type_data = 'Vortices'
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
    resolution = 64
    % resolution = 128
    %resolution = 256;
    % resolution = 512;
    %resolution = 1024;
    % resolution = 2048;
    
    % Resolution of the reference
    resolution_HR = 512;
    % resolution_HR = 1024;
    
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
forcing_type = 'Spring'

% Amplitude of the forcing
ampli_forcing = 10;
% ampli_forcing = 1;

% Frequency of the forcing
freq_f = [3 2]
% freq_f = [0 1];


%% Deterministic subgrid tensor

if nargin == 0
    % Viscosity
    Lap_visco.bool = false;
    
    % % Smagorinsky-like viscosity
    % Smag.bool = false;
    % % HV.bool = false;
    
    % Hyper-viscosity
    HV.bool = true;
    
    if HV.bool
        % HV.order=4;
        HV.order=8;
    end
    
    % Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
    Smag.bool = false;
    
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if Smag.bool
        if Lap_visco.bool
            
            % Use a spatial derivation scheme for the herogeneous
            % disspation
            Smag.spatial_scheme = false;
            
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
            Smag.kappamax_on_kappad = 0.5;
            %Smag.kappamax_on_kappad = 1; % Stable mais petit artefact
            %  d'aliasing  % (better(?))
            
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
plot_moments = false;

% Choose to plot the dissipation by scale
plot_epsilon_k = true;
if sigma.sto & sigma.hetero_energy_flux
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

if nargin == 0
    % Variance tensor a_H
    if stochastic_simulation
        %         if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
        %             sigma.k_c = 0;
        %         else
        %             switch dynamics
        %                 case 'SQG'
        %                     sigma.k_c = 0; % 1/(300 meters)
        %                     % sigma.k_c = 1/(3e2); % 1/(300 meters)
        %                 case '2D'
        %                     error(...
        %                         'The turbulence 2D is not stable under the action of noise');
        %                     %             k_c = 1/(eps); % 1/(100 meters)
        %                 otherwise
        %                     error('Unknown type of dynamics');
        %             end
        %             % a_H is calculated latter in the code using
        %             % a_H = 2 * f_0 / k_c^2
        %             % where f_0 is the Corilis frequency
        %         end
        if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
            if sigma.Smag.bool | ...
                    Lap_visco.bool | ( HV.bool & (HV.order<=4) )
                sigma.kappamin_on_kappamax = 1/2;
                % sigma.kappamin_on_kappamax = 1/4;
                % sigma.kappamin_on_kappamax = 1/8;
            elseif ( HV.bool & (HV.order==8) )
                switch resolution
                    case  128
                        sigma.kappamin_on_kappamax = 1/2;
                    case 64
                        pre=1e-2;
                        sigma.kappamin_on_kappamax = ...
                            (log(1-pre)/log(pre))^(2/HV.order)
                        pre_estim_slope=1e-1;
                        sigma.kappamin_on_kappamax_estim_slope = ...
                            (log(1-pre_estim_slope)/log(pre_estim_slope))...
                            ^(2/HV.order)
                        %                 sigma.kappamin_on_kappamax = 0.45;
                        %                 % sigma.kappamin_on_kappamax = 1/3;
                    otherwise
                        error('unknown');
                end
            else
                warning('kappamin_on_kappamax may be inapropriate');
                sigma.kappamin_on_kappamax = 1/2;
                % sigma.kappamin_on_kappamax = 1/4;
                % sigma.kappamin_on_kappamax = 1/8;
            end
            
            sigma.kappaLS_on_kappamax = 1/8;
        else
            switch resolution
                case  128
                    sigma.kappamin_on_kappamax = 1/2;
                case 64
                    sigma.kappamin_on_kappamax = 1/3;
                otherwise
                    error('unknown');
            end
            
            %             %kappamin_on_kappamax = 1/32;
            %             sigma.kappamin_on_kappamax = 1/2;
            %             % sigma.kappamin_on_kappamax = 1/128;
            %             %         sigma.slope_sigma = - 5;
            %             % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
            
            sigma.kappaLS_on_kappamax = 1/8;
        end
    else
        % If the simulation is deterministic, a_H = 0 and only one simulation
        % is performed
        sigma.k_c = inf; % And then a_H = 0
        N_ech=1;
        plot_moments = false;
    end
    
    % Spectrum slope of sigma dBt
    if sigma.sto
        switch dynamics
            case 'SQG'
                sigma.slope_sigma = - 5/3;
            case '2D'
                sigma.slope_sigma = - 3;
            otherwise
                error('Unknown type of dynamics');
        end
    end
    if  sigma.sto & strcmp(sigma.type_spectrum,'BB')
        sigma.slope_sigma = 0;
        % elseif strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        %     sigma.slope_sigma = nan;
    end
    
    % Rate between the smallest and the largest wave number of sigma dBt
    if sigma.sto
        if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
            sigma.kappamin_on_kappamax = 1/2;
            % sigma.kappamin_on_kappamax = 1/4;
            % sigma.kappamin_on_kappamax = 1/8;
            
            sigma.kappaLS_on_kappamax = 1/8;
        else
            %kappamin_on_kappamax = 1/32;
            sigma.kappamin_on_kappamax = 1/2;
            % sigma.kappamin_on_kappamax = 1/128;
            %         sigma.slope_sigma = - 5;
            % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
            
            sigma.kappaLS_on_kappamax = 1/8;
        end
        
        % Rate between the largest wave number of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappamax_on_kappaShanon = 1;
    end
end

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
if sigma.sto
    eval(['model.sigma.fct_tr_a = @(m,k1,k2) fct_norm_tr_a_theo_' ...
        model.sigma.type_spectrum '(m,k1,k2);']);
end
% eval(['model.sigma.fct_tr_a = @(m,k1,k2,alpha) fct_norm_tr_a_theo_' ...
%     model.sigma.type_spectrum '(m,k1,k2,alpha);']);
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
[~,model_HR] = fct_buoyancy_init(model,resolution_HR);
[~,model] = fct_buoyancy_init(model,resolution);

%% Set up
% if nargin < 1
%     init;
%     day_choose='70';
% end
% if nargin < 2
%     type_data = 'Spectrum'; % 'klein' 'spot' 'Spectrum'
% end
% taille_police = 12;
% HR = true;
% cropped = false;
%
% slop_b_ini=-5/3;
% choice_a = 'choice of k_c'; %  'spectrum_geo' 'averaged_lyapunov' 'kriging' 'choice of k_c';
% eq= 'usual_2D'; % 'transport'  'SQG modif' 'SQG modif deep' 'SQG_modif_ML'
% eq2 = 'bottom ML'; % 'phi tilde C0' 'w C0' 'bottom ML'  'average ML'
% divergent_w = false;
% multiplicative_noise = true;
% plot_each_dt = false;
%  k_c_init = inf;
% % k_c_init = 1/(1e3);
% %k_c_init = 1/(3e2);
% use_save = false;
% day_save = 60;
% advection_duration=3600*24*30*2; % 12 months
plot_modes = true;
plot_error = false;
plot_high_order_moments = false;
nb_modes = 200;
%
% param_SQG.bool=false;
%
% type_adv = 'UQ2';
% % type_adv = 'geostrophic_adv_lagrangian_alone' ; % 'geostrophic_adv_filtered';
% % 'UQ' 'UQ2'
%
% subsamp_x=4;
% subsamp_y=4;
%
% % Colormap
% load('BuYlRd.mat');
% map = BuYlRd; clear BuYlRd
%
% % Version of matlab
% vers = version;
% year = str2double(vers(end-5:end-2));
% subvers = vers(end-1);
% colormap_freeze = ...
%     (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
%
% diff_in_forward=false;
% recompute_a_w=false;
% recompute_kriging_ssh=false;
% recompute_kriging_sst=true;
% mirror=strcmp(type_data,'klein');
% if ~multiplicative_noise
%     warning('no noises');
% end
%
% % N_ech=3;
% N_ech=200;
% bool_HV=true;
%
%
% % Set up parameters
% model.advection.HV.bool=bool_HV;
% model.advection.SQG=param_SQG;
% model.type_data=type_data;
% model.advection.choice_a=choice_a;
% model.advection.meth_anti_alias= 'deriv_LowPass';  % old
% model.advection.N_ech=N_ech;
% model.advection.evol_variance=true;
% model.advection.eq=eq;
% model.advection.eq2=eq2;
% model.advection.multiplicative_noise=multiplicative_noise;
% model.advection.plot_each_dt=plot_each_dt;
% model.advection.k_c_init = k_c_init;
% model.advection.use_save = use_save;
% model.advection.day_save = day_save;
% model.advection.advection_duration=advection_duration;
% model.advection.divergent_w=divergent_w;
% model.advection.nb_modes = nb_modes;
% model.advection.plot_modes = plot_modes;
% model.advection.plot_error = plot_error;



%% Folder with reference
clear subgrid_details
% if model.advection.HV.bool
% add_subgrid_deter = '_HV';
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
% if isinf(model.sigma.k_c) % Deterministic case
model_HR.folder.folder_simu = [ 'images/usual_' model.dynamics ...
    add_subgrid_deter '/' model.type_data ];
% else % Stochastic case
%     model.folder.folder_simu = [ 'images/' model.dynamics ...
%         '_MU' add_subgrid_deter '/' model.type_data ];
% end

if model.advection.forcing.bool
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_forced_turb_' model_HR.advection.forcing.forcing_type];
else
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_free_turb' ];
end
model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/' num2str(model_HR.grid.MX(1)) 'x' num2str(model_HR.grid.MX(2)) ];
model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/Low_Pass_fitlered_version_' ...
    num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];

%% Folder with deterministic model with random initial condition
clear subgrid_details
model_randomIC = model;
model_randomIC.sigma.sto = false;
% if model.advection.HV.bool
% add_subgrid_deter = '_HV';
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
% if isinf(model.sigma.k_c) % Deterministic case
if plot_random_IC
    model_randomIC.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        '_randomIC' add_subgrid_deter '/' model.type_data ];
else
    model_randomIC.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
end
% else % Stochastic case
%     model.folder.folder_simu = [ 'images/' model.dynamics ...
%         '_MU' add_subgrid_deter '/' model.type_data ];
% end
if model.advection.forcing.bool
    model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
        '_forced_turb_' model_randomIC.advection.forcing.forcing_type];
else
    model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
        '_free_turb' ];
end
model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
    '/' num2str(model_randomIC.grid.MX(1)) 'x' num2str(model_randomIC.grid.MX(2)) ];
if plot_random_IC
    if random_IC_large
        model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
            '/large_IC_perturb' ];
    else
        model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
            '/small_IC_perturb' ];
    end
    % else
    %     model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
    %         '/no_IC_perturb' ];
end



%% Folder with deterministic model with random initial condition
clear subgrid_details
model_deter = model;
model_deter.sigma.sto = false;
% if model.advection.HV.bool
% add_subgrid_deter = '_HV';
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
model_deter.folder.folder_simu = [ 'images/usual_' model.dynamics ...
    add_subgrid_deter '/' model.type_data ];
% else % Stochastic case
%     model.folder.folder_simu = [ 'images/' model.dynamics ...
%         '_MU' add_subgrid_deter '/' model.type_data ];
% end
if model.advection.forcing.bool
    model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
        '_forced_turb_' model_deter.advection.forcing.forcing_type];
else
    model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
        '_free_turb' ];
end
model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
    '/' num2str(model_deter.grid.MX(1)) 'x' num2str(model_deter.grid.MX(2)) ];


%% Folder to save plots and files
clear subgrid_details
if model.advection.HV.bool
    add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
elseif model.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
    %add_subgrid_deter = [];
end
if model.sigma.sto & model.sigma.assoc_diff
    add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
end
% if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
%         model.advection.Smag.bool
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool ) | ...
        (model.sigma.sto & model.sigma.Smag.bool )
    % if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
    %  model.advection.Smag.bool
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
    %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
    %         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
    %         '_dealias_ratio_mask_LS_' ...
    %         fct_num2str(model.grid.dealias_ratio_mask_LS)];
    if model.sigma.sto & model.sigma.Smag.bool & ...
            model.sigma.Smag.epsi_without_noise
        add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];
elseif model.sigma.sto & model.sigma.hetero_modulation_V2
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
elseif model.sigma.sto & model.sigma.hetero_energy_flux
    add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
end
if model.sigma.sto & model.sigma.no_noise
    add_subgrid_deter = [add_subgrid_deter '_no_noise'];
end

% if model.sigma.SelfSim_from_LS.bool
%     add_subgrid_deter = [add_subgrid_deter '_SelfSim_from_LS'];
% end

if ~ model.sigma.sto % Deterministic case
    model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
else % Stochastic case
    %     model.folder.folder_simu = [ 'images/' model.dynamics ...
    %         '_MU' add_subgrid_deter '/' model.type_data ];
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_MU' add_subgrid_deter '/' ...
        'type_spectrum_sigma_' model.sigma.type_spectrum '/' ...
        model.type_data ];
end
if model.advection.forcing.bool
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_forced_turb_' model.advection.forcing.forcing_type];
else
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_free_turb' ];
end
model.folder.folder_simu = [ model.folder.folder_simu ...
    '/' num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool)
    subgrid_details = ['kappamax_on_kappad_' ...
        fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
        '_dealias_ratio_mask_LS_' ...
        fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
    if model.advection.Smag.spatial_scheme
        subgrid_details = [ subgrid_details '_spatial_scheme'];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end
if model.sigma.sto
    if model.sigma.Smag.bool
        subgrid_details = ['kappamax_on_kappad_' ...
            fct_num2str(model.sigma.Smag.kappamax_on_kappad) ...
            '_dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if model.sigma.Smag.SS_vel_homo
            subgrid_details = [ subgrid_details '_SS_vel_homo'];
        elseif  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
        if model.advection.Smag.spatial_scheme
            subgrid_details = [ subgrid_details '_spatial_scheme'];
        end
    elseif ( model.sigma.hetero_modulation |  model.sigma.hetero_modulation_V2)
        subgrid_details = ['dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
    end
    %     model.folder.folder_simu = [ model.folder.folder_simu ...
    %         '_kappamin_on_kappamax_' ....
    %         fct_num2str(model.sigma.kappamin_on_kappamax) ];
    if ~ ( exist('subgrid_details','var')==1)
        subgrid_details = [];
    end
    if ~ strcmp(model.sigma.type_spectrum,'EOF')
        subgrid_details = [ subgrid_details ...
            '_kappamin_on_kappamax_' ....
            fct_num2str(model.sigma.kappamin_on_kappamax) ];
        if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
            subgrid_details = [ subgrid_details ...
                '_on_kc_' ....
                fct_num2str(1/model.sigma.k_c) ];
        end
    else
        subgrid_details = [ subgrid_details ...
            '_nbDayLearn_' ...
            fct_num2str(model.sigma.nbDayLearn) ...
            '_Delta_T_on_Delta_t_' ...
            fct_num2str(model.sigma.Delta_T_on_Delta_t) ...;
            '_nb_EOF_' ...
            fct_num2str(model.sigma.nb_EOF)];
    end
    %     if ~ strcmp(model.sigma.type_spectrum,'EOF')
    %         subgrid_details = [ subgrid_details ...
    %             '_kappamin_on_kappamax_' ....
    %             fct_num2str(model.sigma.kappamin_on_kappamax) ];
    %         if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
    %             subgrid_details = [ subgrid_details ...
    %                 '_on_kc_' ....
    %                 fct_num2str(1/model.sigma.k_c) ];
    %         end
    %     end
    %     %     subgrid_details = [ subgrid_details ...
    %     %         '_kappamin_on_kappamax_' ....
    %     %         fct_num2str(model.sigma.kappamin_on_kappamax) ];
    if model.advection.N_ech>1
        subgrid_details = [ subgrid_details ...
            '_N_ech_' ....
            fct_num2str(model.advection.N_ech) ];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end

% Create the folders
fct_create_folder_plots(model,random_IC_large,plot_random_IC)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Comparaison of LU parametrisation
clear subgrid_details
if strcmp(model.sigma.type_spectrum , 'EOF')
    model_SelfSim = model;
    model_SelfSim.sigma.type_spectrum = 'SelfSim_from_LS';
    
    if model_SelfSim.advection.HV.bool
        add_subgrid_deter = ['_HV' '_' fct_num2str(model_SelfSim.advection.HV.order/2)];
    elseif model_SelfSim.advection.Lap_visco.bool
        add_subgrid_deter = '_Lap_visco';
    else
        add_subgrid_deter = '_no_deter_subgrid';
        %add_subgrid_deter = [];
    end
    if model_SelfSim.sigma.sto & model_SelfSim.sigma.assoc_diff
        add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
    end
    % if ( model_SelfSim.advection.HV.bool || model_SelfSim.advection.Lap_visco.bool) && ...
    %         model_SelfSim.advection.Smag.bool
    if ( ( model_SelfSim.advection.HV.bool | model_SelfSim.advection.Lap_visco.bool) & ...
            model_SelfSim.advection.Smag.bool ) | ...
            (model_SelfSim.sigma.sto & model_SelfSim.sigma.Smag.bool )
        % if ( model_SelfSim.advection.HV.bool | model_SelfSim.advection.Lap_visco.bool) & ...
        %  model_SelfSim.advection.Smag.bool
        add_subgrid_deter = [add_subgrid_deter '_Smag'];
        %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
        %         fct_num2str(model_SelfSim.advection.Smag.kappamax_on_kappad) ...
        %         '_dealias_ratio_mask_LS_' ...
        %         fct_num2str(model_SelfSim.grid.dealias_ratio_mask_LS)];
        if model_SelfSim.sigma.sto & model_SelfSim.sigma.Smag.bool & ...
                model_SelfSim.sigma.Smag.epsi_without_noise
            add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
        end
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_modulation
        add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_modulation_V2
        add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_energy_flux
        add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
    end
    if model_SelfSim.sigma.sto & model_SelfSim.sigma.no_noise
        add_subgrid_deter = [add_subgrid_deter '_no_noise'];
    end
    if ~ model_SelfSim.sigma.sto % Deterministic case
        model_SelfSim.folder.folder_simu = [ 'images/usual_' model_SelfSim.dynamics ...
            add_subgrid_deter '/' model_SelfSim.type_data ];
    else % Stochastic case
        %     model_SelfSim.folder.folder_simu = [ 'images/' model_SelfSim.dynamics ...
        %         '_MU' add_subgrid_deter '/' model_SelfSim.type_data ];
        model_SelfSim.folder.folder_simu = [ 'images/' model_SelfSim.dynamics ...
            '_MU' add_subgrid_deter '/' ...
            'type_spectrum_sigma_' model_SelfSim.sigma.type_spectrum '/' ...
            model_SelfSim.type_data ];
    end
    if model_SelfSim.advection.forcing.bool
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '_forced_turb_' model_SelfSim.advection.forcing.forcing_type];
    else
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '_free_turb' ];
    end
    model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
        '/' num2str(model_SelfSim.grid.MX(1)) 'x' num2str(model_SelfSim.grid.MX(2)) ];
    if ( ( model_SelfSim.advection.HV.bool | model_SelfSim.advection.Lap_visco.bool) & ...
            model_SelfSim.advection.Smag.bool)
        subgrid_details = ['kappamax_on_kappad_' ...
            fct_num2str(model_SelfSim.advection.Smag.kappamax_on_kappad) ...
            '_dealias_ratio_mask_LS_' ...
            fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
        if model_SelfSim.advection.Smag.spatial_scheme
            subgrid_details = [ subgrid_details '_spatial_scheme'];
        end
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '/' subgrid_details ];
    end
    if model_SelfSim.sigma.sto
        if model_SelfSim.sigma.Smag.bool
            subgrid_details = ['kappamax_on_kappad_' ...
                fct_num2str(model_SelfSim.sigma.Smag.kappamax_on_kappad) ...
                '_dealias_ratio_mask_LS_' ...
                fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
            if model_SelfSim.sigma.Smag.SS_vel_homo
                subgrid_details = [ subgrid_details '_SS_vel_homo'];
            elseif  model_SelfSim.sigma.proj_free_div
                subgrid_details = [ subgrid_details '_proj_free_div'];
            end
            if model_SelfSim.advection.Smag.spatial_scheme
                subgrid_details = [ subgrid_details '_spatial_scheme'];
            end
        elseif ( model_SelfSim.sigma.hetero_modulation |  model_SelfSim.sigma.hetero_modulation_V2)
            subgrid_details = ['dealias_ratio_mask_LS_' ...
                fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
            if  model_SelfSim.sigma.proj_free_div
                subgrid_details = [ subgrid_details '_proj_free_div'];
            end
        end
        if ~ ( exist('subgrid_details','var')==1)
            subgrid_details = [];
        end
        if ~ strcmp(model_SelfSim.sigma.type_spectrum,'EOF')
            subgrid_details = [ subgrid_details ...
                '_kappamin_on_kappamax_' ....
                fct_num2str(model_SelfSim.sigma.kappamin_on_kappamax) ];
            if strcmp(model_SelfSim.sigma.type_spectrum,'Band_Pass_w_Slope')
                subgrid_details = [ subgrid_details ...
                    '_on_kc_' ....
                    fct_num2str(1/model_SelfSim.sigma.k_c) ];
            end
        else
            subgrid_details = [ subgrid_details ...
                '_nbDayLearn_' ...
                fct_num2str(model_SelfSim.sigma.nbDayLearn) ...
                '_Delta_T_on_Delta_t_' ...
                fct_num2str(model_SelfSim.sigma.Delta_T_on_Delta_t) ...;
                '_nb_EOF_' ...
                fct_num2str(model_SelfSim.sigma.nb_EOF)];
        end
        if model_SelfSim.advection.N_ech>1
            subgrid_details = [ subgrid_details ...
                '_N_ech_' ....
                fct_num2str(model_SelfSim.advection.N_ech) ];
        end
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '/' subgrid_details ];
    end
end

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

dkxdky = (2*pi)^2 / (prod(model.grid.MX.*model.grid.dX));
dkappa = sqrt(dkxdky);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
N_ech=model.advection.N_ech;

%%

t_last_plot = -inf;
day_last_plot = - inf;
model.advection.step='finite_variation';

% t_ini=1;

folder_ref = model.folder.folder_simu;
name_file = [model.folder.folder_simu '/files/' num2str(first_day) '.mat'];
load(name_file)
model.folder.folder_simu = folder_ref;

dt=model.advection.dt_adv;
N_t = ceil(model.advection.advection_duration/dt);

% x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% % if model.mirror
% %     y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
% % else
% %     y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% % end
% My=model.grid.MX(2);
% % if model.mirror
% %     My=My/2;
% % end

F_save = [];
F_save2 = [];
bt1_HR_vect = [];
bt1_LR_vect = [];

error_vs_t = [];
error_vs_t_SelfSim = [];
 v_day = [];

t_ini = first_day*24*3600/dt;
trigger = false;
dt_loop = dt;
%t_ini=1700000
for t_loop=t_ini:N_t
    %     %% Plot
    % t_loop=1;
    day_num = (floor(t_loop*dt_loop/24/3600));
    
    if day_num > day_last_plot
        % if (t_loop - t_last_plot)*dt >= 3600*24*1
        day_num = (floor(t_loop*dt_loop/24/3600));
        day = num2str(day_num);
        day
        day_last_plot = day_num;
        
        % if ~(exist('time','var')==1)
        time =t_loop*dt_loop;
        % end
        
        model.advection.plot_modes = plot_modes;
        model.advection.nb_modes = nb_modes;
        t_last_plot=t_loop;
        id_part=1;
        
        width=1.2e3;
        height=0.5e3;
        if strcmp(model.type_data,'klein')
            width=width/2;
            r_c_ax = 0.5;
        else
            r_c_ax =1/1.5;
            %             r_c_ax =1;
        end
        %         X0=0.5e3*[1 1];
        X0=[0 1];
        
        %% Specific day
        %         warning('the day is specified manually');
        %         day = day_choose;
        %         day = num2str(day);
        
        %% Load meth with random IC
        clear fft_b;
        if (model.advection.N_ech>1)
            %         if plot_random_IC & (model.advection.N_ech>1)
            name_file_randomIC = [model_randomIC.folder.folder_simu ...
                '/files/' num2str(day) '.mat'];
            load(name_file_randomIC,'fft_T_adv_part','fft_buoy_part','fft_b');
            if (exist('fft_b','var')==1)
            elseif (exist('fft_buoy_part','var')==1)
                fft_b = fft_buoy_part;
            elseif (exist('fft_T_adv_part','var')==1)
                fft_b = fft_T_adv_part;
            else
                error('Cannot find buoyancy field')
            end
            fft_b_classic = fft_b;clear fft_b;
        end
        
        %% Load determinsitic meth without random IC
        clear fft_b;
        %         if (model.advection.N_ech>1)
        name_file_deter = [model_deter.folder.folder_simu ...
            '/files/' num2str(day) '.mat'];
        load(name_file_deter,'fft_T_adv_part','fft_buoy_part','fft_b');
        if (exist('fft_b','var')==1)
        elseif (exist('fft_buoy_part','var')==1)
            fft_b = fft_buoy_part;
        elseif (exist('fft_T_adv_part','var')==1)
            fft_b = fft_T_adv_part;
        else
            error('Cannot find buoyancy field')
        end
        fft_b_deter = fft_b;clear fft_b;
        model_deter.grid.x=model_deter.grid.x_ref;
        model_deter.grid.y=model_deter.grid.y_ref;
        model_deter.folder.colormap=model.folder.colormap;
        %         end
        
        %% Load selfSim
        if strcmp(model.sigma.type_spectrum , 'EOF')
            model_SelfSim_ref = model_SelfSim;
            name_file_SelfSim = ...
                [model_SelfSim.folder.folder_simu '/files/' num2str(day) '.mat'];
            clear fft_b_SelfSim fft_b fft_buoy_part;
            if exist( name_file_SelfSim,'file')==2
                % clear fft_b fft_buoy_part;
                load(name_file_SelfSim,'fft_T_adv_part','fft_buoy_part','fft_b');
                % model_SelfSim.folder.folder_simu = folder_SelfSim_ref;
                if ~(exist('fft_b','var')==1)
                    fft_b_SelfSim = fft_buoy_part; clear fft_buoy_part
                else
                    fft_b_SelfSim = fft_b; clear fft_b
                end
            else
                warning('Cannot find the following file');
                fprintf([ name_file ' \n']);
            end
            
            
        end
        
        %% Load
        model_ref = model;
        name_file = [model.folder.folder_simu '/files/' num2str(day) '.mat'];
        clear fft_b fft_buoy_part;
        if exist( name_file,'file')==2
            % clear fft_b fft_buoy_part;
            load(name_file);
            model.folder.folder_simu = folder_ref;
            if ~(exist('fft_b','var')==1)
                fft_b = fft_buoy_part;
            end
        else
            warning('Cannot find the following file');
            fprintf([ name_file ' \n']);
            return
        end
        
        %model.odg_b = model.odg_b*3;
        %
        %         fct_plot(model,fft_buoy_part,day)
        
        %% Load reference
        name_file_HR = [model_HR.folder.folder_simu '/files/' num2str(day) '.mat'];
        load(name_file_HR,'fft_buoy_part_ref','spectrum_ref');
        model_HR.grid.x=model_HR.grid.x_ref;
        model_HR.grid.y=model_HR.grid.y_ref;
        model_HR.folder.colormap=model.folder.colormap;
        
        %%
        
        fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
        
        %%
        %         if model.advection.Smag.bool
        %             [coef_diff_aa,coef_diff] = fct_coef_diff(model,fft_buoy_part);
        %             figure(9);
        %             subplot(1,2,1);
        %             imagesc(model.grid.x_ref,model.grid.y_ref,...
        %                 real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
        %             subplot(1,2,2);
        %             imagesc(model.grid.x_ref,model.grid.y_ref,...
        %                 coef_diff_aa');axis xy;axis equal;colorbar;
        %             drawnow
        %             eval( ['print -depsc ' model.folder.folder_simu ...
        %                 '/dissip_coef/' day '.eps']);
        %         end
        %%
        % Plots
        
        % Colormap
        load('BuYlRd.mat');
        model.folder.colormap = BuYlRd; clear BuYlRd
        
        % Version of matlab
        vers = version;
        year = str2double(vers(end-5:end-2));
        subvers = vers(end-1);
        model.folder.colormap_freeze = ...
            (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
        %%
        
        if model.advection.N_ech > 1
            model.advection.plot_moments = true;
            fct_plot_post_process(model,fft_b,day);
            pause(0.1);
        end
        
        fct_plot_post_process(model_deter,fft_b_deter,day);
        %         fct_plot_post_process(model_HR,fft_buoy_part_ref,day);
        
        model.advection.plot_moments = false;
        [spectrum,name_plot] = fct_plot_post_process(model,fft_b,day);
        if strcmp(model.sigma.type_spectrum , 'EOF')
            color_SelfSim = [1 0 1];
            fct_spectrum_multi(model,fft_b_SelfSim,color_SelfSim);
            fct_spectrum_multi(model,fft_b_deter,'m');
            fct_spectrum_multi(model,fft_buoy_part_ref,'r');
            legend('-5/3',...
                ['LU EOF' num2str(resolution)],...
                ['LU Self.Sim.' num2str(resolution)],...
                ['Deter. ' num2str(resolution)],...
                ['Deter. ' num2str(resolution_HR')]);
        else
            fct_spectrum_multi(model,fft_b_deter,'m');
            fct_spectrum_multi(model,fft_buoy_part_ref,'r');
            legend('-5/3',['LU ' num2str(resolution)],...
                ['Deter. ' num2str(resolution)],...
                ['Deter. ' num2str(resolution_HR')]);
            %legend('-5/3','LU 128','Deter. 128','Deter. 1024')
        end
        eval( ['print -depsc ' model.folder.folder_simu '/Spectrum/' day '.eps']);
        %         if model.advection.plot_dissip
        %             fct_plot_dissipation(model,fft_buoy_part,sigma_on_sq_dt,day);
        %         end
        
        %% Distance to reference
        % Spectrum discrepancy
        temp(1,1) = ...
            dkappa * sum(abs(spectrum(:) - spectrum_ref(:))) ;
        % Error
        error2 = 1/prod(model.grid.MX)^2 * ...
            sum(sum(abs(bsxfun(@plus,fft_b, ...
            - fft_buoy_part_ref)).^2 , 1) ,2);
        % Bias
        bias2 = 1/prod(model.grid.MX)^2 * ...
            sum(sum(abs(bsxfun(@plus,mean(fft_b ,4),...
            - fft_buoy_part_ref)).^2 , 1) ,2);
        temp(1,2) = sqrt(bias2);
        % RMSE
        temp(1,3) = sqrt(mean(error2 ,4));
        % min distance
        temp(1,4) = min(sqrt(error2),[],4);
        % Concatenation
        error_vs_t = [ error_vs_t ; temp ];
        v_day = [ v_day eval(day)];
        
        figure1111=figure(1111);
        close(figure1111)
        figure1111=figure(1111);
        if strcmp(model.sigma.type_spectrum , 'EOF')
            % Spectrum discrepancy
            temp(1,1) = nan ;
            % Error
            error2 = 1/prod(model.grid.MX)^2 * ...
                sum(sum(abs(bsxfun(@plus,fft_b_SelfSim, ...
                - fft_buoy_part_ref)).^2 , 1) ,2);
            % Bias
            bias2 = 1/prod(model.grid.MX)^2 * ...
                sum(sum(abs(bsxfun(@plus,mean(fft_b_SelfSim ,4),...
                - fft_buoy_part_ref)).^2 , 1) ,2);
            temp(1,2) = sqrt(bias2);
            % RMSE
            temp(1,3) = sqrt(mean(error2 ,4));
            % min distance
            temp(1,4) = min(sqrt(error2),[],4);
            % Concatenation
            error_vs_t_SelfSim = [ error_vs_t_SelfSim ; temp ];
            hold on;
            %figure;
            plot(v_day ,[error_vs_t(:,2:4) error_vs_t_SelfSim(:,2:4)]');
            %hold off;
            
            legend('Bias EOF','RMSE EOF','Min. dist. EOF',...
                'Bias SelfSim','RMSE SelfSim','Min. dist. SelfSim');
        else
            figure;plot(v_day ,error_vs_t(:,2:4)');
            legend('Bias','RMSE','Min. dist.');
        end
        drawnow;
        eval( ['print -depsc ' model.folder.folder_simu ...
            '/error_along_time.eps']);
        
        if model.advection.N_ech > 1
            if strcmp(model.sigma.type_spectrum , 'EOF')
                plot_error_ensemble_comp_EOF_SelfSim
            else
                plot_error_ensemble
            end
        end
        
        %%
        if strcmp( model.sigma.type_spectrum, 'SelfSim_from_LS')
            fft_w = SQG_large_UQ(model, fft_b(:,:,:,1));
            fct_sigma_spectrum_abs_diff_postprocess(...
                model,fft_w,true,day);
        end
        
    end
end

% Square root
% error_vs_t = sqrt( error_vs_t);

end