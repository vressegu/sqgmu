function plot_post_process_2(stochastic_simulation,type_data,resolution,forcing, ...
    sigma,Lap_visco,HV,Smag)
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


if nargin == 0
    
    % Deterministic or random model
    stochastic_simulation = false;
    sigma.sto = stochastic_simulation;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
    
    if sigma.sto
        % Type of spectrum for sigma dBt
        % type_spectrum = 'Band_Pass_w_Slope'; % as in GAFD part II
        %type_spectrum = 'Low_Pass_w_Slope';
        % Spectrum cst for k<km ans slope for k>km
        % type_spectrum = 'Low_Pass_streamFct_w_Slope';
        % Matern covariance for the streamfunction
        % spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
        % ~ k2 for k<km ans slope for k>km
        % type_spectrum = 'BB';
        % type_spectrum = 'Bidouille';
        sigma.type_spectrum = 'SelfSim_from_LS';
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
                | sigma.hetero_modulation_V2
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
            
            % Use a spatial derivation scheme for the herogeneous
            % disspation
            Smag.spatial_scheme = false;
            
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
            warning('There isno noise here');
        end
    end
end

% Number of realizations in the ensemble
if nargin == 0
    N_ech=200;
else
    N_ech=1;
end
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )

% Duration of the simulation (in seconds)
advection_duration = 3600*24*30;
%advection_duration = 3600*24*1000;
% % advection_duration = 3600*24*20; % 20 days

if nargin == 0
    % Type of initial condtions
    type_data = 'Vortices';
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
    % resolution = 128;
    %resolution = 256;
    % resolution = 512;
    resolution = 1024;
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
    % end
    
    % Type de forcing
    % forcing_type = 'Kolmogorov';
    forcing_type = 'Spring';
    
    % Amplitude of the forcing
    ampli_forcing = 10;
    % ampli_forcing = 1;
    
    % Frequency of the forcing
    freq_f = [3 2];
    % freq_f = [0 1];
end

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
            % Smag.dealias_ratio_mask_LS = 1/8;
            %     Smag.dealias_ratio_mask_LS = 1/8;
            %     %Smag.dealias_ratio_mask_LS = 1/2;
            Smag.dealias_ratio_mask_LS = 1;
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
plot_epsilon_k = false;
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
        if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
            sigma.k_c = 0;
        else
            switch dynamics
                case 'SQG'
                    sigma.k_c = 1/(3e2); % 1/(300 meters)
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
        end
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
    if  sigma.sto & strcmp(sigma.type_spectrum,'BB')
        sigma.slope_sigma = 0;
        % elseif strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        %     sigma.slope_sigma = nan;
    end
    
    % <<<<<<< HEAD
    % % Rate between the smallest and the largest wave number of sigma dBt
    % if sigma.sto & strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
    %     sigma.kappamin_on_kappamax = 1/2;
    %     % sigma.kappamin_on_kappamax = 1/4;
    %     % sigma.kappamin_on_kappamax = 1/8;
    %
    %     sigma.kappaLS_on_kappamax = 1/8;
    % else
    %     %kappamin_on_kappamax = 1/32;
    %     sigma.kappamin_on_kappamax = 1/2;
    %     % sigma.kappamin_on_kappamax = 1/128;
    %     %         sigma.slope_sigma = - 5;
    %     % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
    % =======
    % if sigma.sto
    %     % Rate between the smallest and the largest wave number of sigma dBt
    %     if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
    %         sigma.kappamin_on_kappamax = 1/2;
    %         % sigma.kappamin_on_kappamax = 1/4;
    %         % sigma.kappamin_on_kappamax = 1/8;
    %
    %         sigma.kappaLS_on_kappamax = 1/8;
    %     else
    %         %kappamin_on_kappamax = 1/32;
    %         sigma.kappamin_on_kappamax = 1/2;
    %         % sigma.kappamin_on_kappamax = 1/128;
    %         %         sigma.slope_sigma = - 5;
    %         % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
    %
    %         sigma.kappaLS_on_kappamax = 1/8;
    %     end
    %     % >>>>>>> 66e0d274c3dc3d35c2b1138c12b022fd3a0806f5
    %
    %     % Rate between the largest wave number of sigma dBt and the largest wave
    %     % number of the simulation
    %     sigma.kappamax_on_kappaShanon = 1;
    % end
    % if nargin == 0
    % Rate between the smallest and the largest wave number of sigma dBt
    if sigma.sto
        switch sigma.type_spectrum
            case {'SelfSim_from_LS','EOF','Euler_EOF'}
                
                pre_estim_slope=1e-1;
                %%
                pre_5 = 5e-2;
                sigma.kappamin_on_kappamax = ...
                    (log(1-pre_5)/log(pre_estim_slope))^(2/HV.order);
                %         sigma.kappamin_on_kappamax = ...
                %             (log(1-pre_estim_slope)/log(pre_estim_slope))^(2/HV.order);
                %%
                %                 pre=1e-2;
                %                 sigma.kappamin_on_kappamax = ...
                %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order);
                %%
                sigma.kappamin_on_kappamax_estim_slope = ...
                    (log(1-pre_estim_slope)/log(pre_estim_slope))...
                    ^(2/HV.order);
                
                sigma.kappaLS_on_kappamax = 1/8;
                
            otherwise
                switch resolution
                    case  128
                        sigma.kappamin_on_kappamax = 1/2;
                    case 64
                        sigma.kappamin_on_kappamax = 1/3;
                    otherwise
                        error('unknown');
                end
        end
        
        %         if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
        %             if sigma.Smag.bool | ...
        %                     Lap_visco.bool | ( HV.bool & (HV.order<=4) )
        %                 sigma.kappamin_on_kappamax = 1/2;
        %                 % sigma.kappamin_on_kappamax = 1/4;
        %                 % sigma.kappamin_on_kappamax = 1/8;
        %             elseif ( HV.bool & (HV.order==8) )
        %                 switch resolution
        %                     case  128
        %                 % sigma.kappamin_on_kappamax = 1/2;
        %                 pre=1e-2;
        %                 pre_estim_slope=1e-1;
        %                 sigma.kappamin_on_kappamax = ...
        %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order)
        %                 sigma.kappamin_on_kappamax_estim_slope = ...
        %                     (log(1-pre_estim_slope)/log(pre_estim_slope))...
        %                     ^(2/HV.order)
        %                     case 64
        %                 pre=1e-2;
        %                 pre_estim_slope=1e-1;
        %                 sigma.kappamin_on_kappamax = ...
        %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order)
        %                 sigma.kappamin_on_kappamax_estim_slope = ...
        %                     (log(1-pre_estim_slope)/log(pre_estim_slope))...
        %                     ^(2/HV.order)
        % %                 sigma.kappamin_on_kappamax = 0.45;
        % %                 % sigma.kappamin_on_kappamax = 1/3;
        %                     otherwise
        %                         error('unknown');
        %                 end
        %             else
        %                 warning('kappamin_on_kappamax may be inapropriate');
        %                 sigma.kappamin_on_kappamax = 1/2;
        %                 % sigma.kappamin_on_kappamax = 1/4;
        %                 % sigma.kappamin_on_kappamax = 1/8;
        %             end
        %
        %             sigma.kappaLS_on_kappamax = 1/8;
        %         else
        %             switch resolution
        %                 case  128
        %                     sigma.kappamin_on_kappamax = 1/2;
        %                 case 64
        %                     sigma.kappamin_on_kappamax = 1/3;
        %                 otherwise
        %                     error('unknown');
        %             end
        %
        % %             %kappamin_on_kappamax = 1/32;
        % %             sigma.kappamin_on_kappamax = 1/2;
        % %             % sigma.kappamin_on_kappamax = 1/128;
        % %             %         sigma.slope_sigma = - 5;
        % %             % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
        %
        %             sigma.kappaLS_on_kappamax = 1/8;
        %         end
        
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
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

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

%% Folder to save plots and files
% if model.advection.HV.bool
%     add_subgrid_deter = '_HV';
% elseif model.advection.Lap_visco.bool
%     add_subgrid_deter = '_Lap_visco';
% else
%     add_subgrid_deter = '_no_deter_subgrid';
%     %add_subgrid_deter = [];
% end
% if model.sigma.sto & model.sigma.assoc_diff
%     add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
% end
% % if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
% %         model.advection.Smag.bool
% if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%         model.advection.Smag.bool ) | model.sigma.Smag.bool
%     % if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%     %  model.advection.Smag.bool
%     add_subgrid_deter = [add_subgrid_deter '_Smag'];
%     %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
%     %         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
%     %         '_dealias_ratio_mask_LS_' ...
%     %         fct_num2str(model.grid.dealias_ratio_mask_LS)];
%     if model.sigma.sto & model.sigma.Smag.bool & ...
%             model.sigma.Smag.epsi_without_noise
%         add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
%     end
% elseif model.sigma.sto & model.sigma.hetero_modulation
%     add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];
% elseif model.sigma.sto & model.sigma.hetero_modulation_V2
%     add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
% elseif model.sigma.sto & model.sigma.hetero_energy_flux
%     add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
% end
% if model.sigma.no_noise
%     add_subgrid_deter = [add_subgrid_deter '_no_noise'];
% end
%
% % if model.sigma.SelfSim_from_LS.bool
% %     add_subgrid_deter = [add_subgrid_deter '_SelfSim_from_LS'];
% % end
%
% if ~ model.sigma.sto % Deterministic case
%     model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
%         add_subgrid_deter '/' model.type_data ];
% else % Stochastic case
%     %     model.folder.folder_simu = [ 'images/' model.dynamics ...
%     %         '_MU' add_subgrid_deter '/' model.type_data ];
%     model.folder.folder_simu = [ 'images/' model.dynamics ...
%         '_MU' add_subgrid_deter '/' ...
%         'type_spectrum_sigma_' model.sigma.type_spectrum '/' ...
%         model.type_data ];
% end
% if model.advection.forcing.bool
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '_forced_turb' ];
% else
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '_free_turb' ];
% end
% model.folder.folder_simu = [ model.folder.folder_simu ...
%     '/' num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];
% if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%         model.advection.Smag.bool)
%     subgrid_details = ['kappamax_on_kappad_' ...
%         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
%         '_dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% end
% if model.sigma.sto & model.sigma.Smag.bool
%     subgrid_details = ['kappamax_on_kappad_' ...
%         fct_num2str(model.sigma.Smag.kappamax_on_kappad) ...
%         '_dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     if model.sigma.Smag.SS_vel_homo
%         subgrid_details = [ subgrid_details '_SS_vel_homo'];
%     elseif  model.sigma.proj_free_div
%         subgrid_details = [ subgrid_details '_proj_free_div'];
%     end
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% elseif model.sigma.sto & ...
%         ( model.sigma.hetero_modulation |  model.sigma.hetero_modulation_V2)
%     subgrid_details = ['dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     if  model.sigma.proj_free_div
%         subgrid_details = [ subgrid_details '_proj_free_div'];
%     end
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% end
% if model.sigma.sto
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '_kappamin_on_kappamax_' ....
%         fct_num2str(model.sigma.kappamin_on_kappamax) ];
%     subgrid_details = [ subgrid_details ...
%         '_kappamin_on_kappamax_' ....
%         fct_num2str(model.sigma.kappamin_on_kappamax) ];
% end
%
% % Create the folders
% fct_create_folder_plots(model)
%
% % Colormap
% load('BuYlRd.mat');
% model.folder.colormap = BuYlRd; clear BuYlRd
%
% % Version of matlab
% vers = version;
% year = str2double(vers(end-5:end-2));
% subvers = vers(end-1);
% model.folder.colormap_freeze = ...
%     (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Folder to save plots and files
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
        '_forced_turb' ];
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
%
% if model.sigma.sto & model.sigma.Smag.bool
%     subgrid_details = ['kappamax_on_kappad_' ...
%         fct_num2str(model.sigma.Smag.kappamax_on_kappad) ...
%         '_dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     if model.sigma.Smag.SS_vel_homo
%         subgrid_details = [ subgrid_details '_SS_vel_homo'];
%     elseif  model.sigma.proj_free_div
%         subgrid_details = [ subgrid_details '_proj_free_div'];
%     end
%     if model.advection.Smag.spatial_scheme
%         subgrid_details = [ subgrid_details '_spatial_scheme'];
%     end
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% elseif model.sigma.sto & ...
%         ( model.sigma.hetero_modulation |  model.sigma.hetero_modulation_V2)
%     subgrid_details = ['dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     if  model.sigma.proj_free_div
%         subgrid_details = [ subgrid_details '_proj_free_div'];
%     end
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% end
% if model.sigma.sto
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '_kappamin_on_kappamax_' ....
%         fct_num2str(model.sigma.kappamin_on_kappamax) ];
%     subgrid_details = [ subgrid_details ...
%         '_kappamin_on_kappamax_' ....
%         fct_num2str(model.sigma.kappamin_on_kappamax) ];
% end
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
    elseif ( model.sigma.hetero_modulation ...
            |  model.sigma.hetero_modulation_V2 ...
            |  model.sigma.hetero_modulation_Smag ...
            | model.sigma.hetero_energy_flux )
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
    subgrid_details = [ subgrid_details ...
        '_kappamin_on_kappamax_' ....
        fct_num2str(model.sigma.kappamin_on_kappamax) ];
    if model.advection.N_ech > 1
        subgrid_details = [ subgrid_details ...
            '_N_ech_' ....
            fct_num2str(model.advection.N_ech) ];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end

% Create the folders
% fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
N_ech=model.advection.N_ech;

%%

t_last_plot = -inf;
model.advection.step='finite_variation';

t_ini=1;

folder_ref = model.folder.folder_simu;
name_file = [model.folder.folder_simu '/files/' num2str(0) '.mat'];
load(name_file)
model.folder.folder_simu = folder_ref;

dt=model.advection.dt_adv;
N_t = ceil(model.advection.advection_duration/dt);

x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% if model.mirror
%     y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
% else
%     y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% end
My=model.grid.MX(2);
% if model.mirror
%     My=My/2;
% end

F_save = [];
F_save2 = [];
bt1_HR_vect = [];
bt1_LR_vect = [];
v_epsilon_dissip = [];
v_epsilon_th_scale = [];
v_time = [];


% trigger = false;
% %t_ini=1700000
% for t_loop=t_ini:N_t
%     %     %% Plot
%     % t_loop=1;
%     if (t_loop - t_last_plot)*dt >= 3600*24*1
%         day = num2str(floor(t_loop*dt/24/3600));
%         t_last_plot = t_loop;
%         model.advection.plot_modes = plot_modes;
%         model.advection.nb_modes = nb_modes;
%         t_last_plot=t_loop;
%         id_part=1;

day_last_plot = -inf;
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
        
        %% Load
        model_ref = model;
        name_file = [model.folder.folder_simu '/files/' num2str(day) '.mat'];
        clear fft_b fft_buoy fft_buoy_part fft_T_adv_part
        if exist( name_file,'file')==2
            load(name_file)
            model.folder.folder_simu = folder_ref;
        else
            warning('Cannot find the following file');
            fprintf([ name_file ' \n']);
            return
        end
        
        %model.odg_b = model.odg_b*3;
        %
        %         fct_plot(model,fft_buoy_part,day)
        %%
        
        %fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        
        %%
        % if ~(exist('time','var')==1)
        if (exist('t','var')==1)
            time =t*dt;
        end
        %         if ~(exist('fft_b','var')==1)
        %             fft_b =fft_buoy_part;
        %         end
        if (exist('fft_b','var')==1)
        elseif (exist('fft_buoy_part','var')==1)
            fft_b = fft_buoy_part;
        elseif (exist('fft_T_adv_part','var')==1)
            fft_b = fft_T_adv_part;
        else
            error('Cannot find buoyancy field')
        end
        
        if ~isfield(model.sigma,'sto')
            model.sigma.sto = (model.sigma.a0>0);
        end
        fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
        a_0_LS = mean(sigma_dBt_on_sq_dt(:).^2);
        %         a_0_LS = mean(sigma_dBt_dt(:).^2)*model.advection.dt_adv;
        %         %a_0_LS = mean(sigma_dBt_dt(:).^2)*model.advection.dt_adv/2;
        if model.sigma.sto
            a_0_LS
        end
        
        %%
        if model.advection.Smag.bool || ...
                ( model.sigma.sto && model.sigma.Smag.bool) ...
                || ( model.sigma.sto && ( ...
                model.sigma.hetero_modulation ...
                || model.sigma.hetero_modulation_V2 ))
            id_part=1;
            % Coefficient coef_Smag to target a specific diffusive scale
            if model.advection.Smag.bool
                [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                    fft_b(:,:,:,id_part));
                coef_diff_aa = ...
                    model.advection.Smag.coef_Smag * coef_diff_aa ;
                coef_diff = ...
                    model.advection.Smag.coef_Smag * coef_diff ;
            elseif model.sigma.Smag.bool
                [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                    fft_b(:,:,:,id_part));
                if model.sigma.a0_SS > eps
                    coef_diff = ...
                        (1 + model.sigma.a0_LS / model.sigma.a0_SS) * ...
                        model.sigma.Smag.coef_Smag * coef_diff ;
                    coef_diff_aa = ...
                        (1 + model.sigma.a0_LS / model.sigma.a0_SS) * ...
                        model.sigma.Smag.coef_Smag * coef_diff_aa ;
                elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                    % The absolute diffusivity diagnosed from the large-scale
                    % kinematic spectrum is too weak. It suggests that there
                    % are few small scales and no subgrid terms is needed.
                    % Moreover, setting subgris terms to zero prevent numerical
                    % errors.
                    coef_diff = 0;
                    coef_diff_aa = 0;
                else
                    error('Unknow case');
                end
                
            elseif model.sigma.hetero_modulation ...
                    | model.sigma.hetero_modulation_V2
                [coef_diff_aa,coef_diff] = ...
                    fct_coef_estim_AbsDiff_heterogeneous(...
                    model,fft_w(:,:,:,id_part));
                coef_diff_aa = ...
                    model.sigma.a0(id_part)/2 * coef_diff_aa ;
                coef_diff = ...
                    model.sigma.a0(id_part)/2 * coef_diff ;
            elseif model.sigma.hetero_energy_flux
                coef_diff_aa = ...
                    fct_epsilon_k_onLine(model,fft_b,fft_w);
                coef_diff_aa = ...
                    model.sigma.a0(id_part)/2 * coef_diff_aa ;
            end
            figure(9);
            subplot(1,2,1);
            imagesc(model.grid.x_ref,model.grid.y_ref,...
                real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
            subplot(1,2,2);
            imagesc(model.grid.x_ref,model.grid.y_ref,...
                coef_diff_aa');axis xy;axis equal;colorbar;
            drawnow
            eval( ['print -depsc ' model.folder.folder_simu ...
                '/dissip_coef/' day '.eps']);
        end
        %    if model.advection.Smag.bool || model.sigma.Smag.bool ...
        %              || model.sigma.hetero_modulation
        %  % Coefficient coef_Smag to target a specific diffusive scale
        %         if model.sigma.Smag.bool
        %  [coef_diff_aa,coef_diff] = fct_coef_diff(model,fft_b);
        %   coef_diff_aa = model.sigma.Smag.coef_Smag * coef_diff_aa ;
        %          else
        %   coef_diff_aa = model.advection.Smag.coef_Smag * coef_diff_aa ;
        %          end
        %         figure(9);
        %        subplot(1,2,1);
        %         imagesc(model.grid.x_ref,model.grid.y_ref,...
        %         real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
        %                 subplot(1,2,2);
        %                 imagesc(model.grid.x_ref,model.grid.y_ref,...
        %                     coef_diff_aa');axis xy;axis equal;colorbar;
        %                 drawnow
        %        eval( ['print -depsc ' model.folder.folder_simu ...
        %                     '/dissip_coef/' day '.eps']);
        %             end
        %%
        % Plots
        [spectrum,name_plot,int_epsilon] = ...
            fct_plot_post_process(model,fft_b,day);
        %    fct_plot(model,fft_b,day);
        
        if model.advection.plot_dissip
            if ~model.sigma.sto
                sigma_on_sq_dt = 0;
            end
            epsilon_dissip = fct_plot_dissipation(model,fft_b,sigma_on_sq_dt,day);
            v_epsilon_dissip = [ v_epsilon_dissip epsilon_dissip];
            v_time = [ v_time time];
        end
        
        if model.sigma.sto & ...
                strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            fct_sigma_spectrum_abs_diff(model,fft2(w),true,day);
            %    sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) ...
            %      * sigma; clear sigma
            %     model.sigma.a0 = a0;
            %    model.sigma.a0_on_dt = ...
            %        model.sigma.a0 / model.advection.dt_adv;
            %     % Diffusion coefficient
            %    model.advection.coef_diff = 1/2 * model.sigma.a0;
            %        warning('The CFL should be changed');
        end
        
        
        if model.sigma.sto & ...
                ( model.sigma.Smag.bool | model.sigma.assoc_diff )
            slope_sigma=model.sigma.slope_sigma
            a0_LS=model.sigma.a0_LS
            a0_SS=model.sigma.a0_SS
        end
        
        %         if model.advection.cov_and_abs_diff
        %             abs_diff = sum(cov_w(t_ref_cov:end))*model.advection.dt_adv;
        %             figure(36)
        %             plot(model.advection.dt_adv*(0:(N_t-1))/3600/24,cov_w);
        %             hold on;
        %             plot(t_ref_cov*[1 1]*model.advection.dt_adv/3600/24,...
        %                 max(abs(cov_w(~isnan(cov_w))))*[-1 1],'r');
        %             hold off
        %         end
        
        % Dissipation by scale
        if plot_epsilon_k
            %if model.advection.plot_epsilon_k
            epsilon_th_scales = fct_plot_epsilon_k(model,fft_b,day);
            % fct_plot_epsilon_k(model,fft_b,int_epsilon,day);
            v_epsilon_th_scale = [ v_epsilon_th_scale epsilon_th_scales];
        end
        dt = model.advection.dt_adv
        
        
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
        %         %%
        %         % Plots
        %         [spectrum,name_plot,int_epsilon] = ...
        %                 fct_plot(model,fft_buoy_part,day);
        %
        %         if model.advection.plot_dissip
        %             fct_plot_dissipation(model,fft_buoy_part,sigma_on_sq_dt,day);
        %         end
    end
    
    %     % Dissipation by scale
    %     if (t_loop == t_last_plot + 1) &  plot_epsilon_k
    %         fct_plot_epsilon_k(model,fft_buoy_part,int_epsilon,day);
    %     end
end


if plot_epsilon_k && model.advection.plot_dissip
    figure(88);plot(v_time,v_epsilon_th_scale,'r');
    hold on;plot(v_time,v_epsilon_dissip,'b');
    hold off
    eval( ['print -depsc ' model.folder.folder_simu '/epsilon_vs_time.eps']);
end


end