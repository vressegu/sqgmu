function post_process_error_grid_ref(stochastic_simulation,...
    type_data,resolution,forcing,Lap_visco,HV,Smag)
    %Lap_visco,HV,Smag,day_choose)
% plot the same thing that fct_fft_advection_sto
%

%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
% init;

%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
% dynamics = '2D';

if nargin == 0
    % Deterministic or random model
    stochastic_simulation = false;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
end

% Duration of the simulation (in seconds)
advection_duration = 3600*24*30;
% advection_duration = 3600*24*1000;
% % advection_duration = 3600*24*20; % 20 days

% Number of realizations in the ensemble
N_ech=1;
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )

if nargin == 0
    % Type of initial condtions
    type_data ='Constantin_case2' ;
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
    %resolution = 128;
    %resolution = 512;
    %resolution = 128;
    resolution = 1024;
    %resolution = 2048;
    
    resolution_LR = 128;
    
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
    
    if Lap_visco.bool
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1/8;
        %     dealias_ratio_mask_LS = 1/8;
        %     %dealias_ratio_mask_LS = 1/2;
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
        %     %     Smag.kappamax_on_kappad = 2;
        %     % % Smag.kappamax_on_kappad = 1.5; % Better, still small oscillations or just pixels?
        %     % %    %  Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
        %     Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
        %     %  d'aliasing
        Smag.kappamax_on_kappad = 1; % Stable mais petit artefact
        %  d'aliasing
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Factor in front of the additional constant dissipation
        % Set to 0 for no additional constant dissipation
        %    HV.weight_cst_dissip = 1/10;
        %     HV.weight_cst_dissip = 1/10; % no aliasing
        Smag.weight_cst_dissip = 0;
    elseif HV.bool
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
        Smag.kappamax_on_kappad = 1.1;% still small oscillations or just pixels?
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Factor in front of the additional constant dissipation
        % Set to 0 for no additional constant dissipation
        % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
        %     % HV.weight_cst_dissip = 1/3; % still a bit of (stable) aliasing
        %     HV.weight_cst_dissip = 1/3;
        %     % HV.weight_cst_dissip = 0;
        
        Smag.weight_cst_dissip = 1/1;
        % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
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
plots_bool = true;

% Choose to plot one-point one-time moments each day
plot_moments = false;

% Plot dissipations terms
plot_dissip = false;

% Begin simulation from a precomputed field?
use_save = true;
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
    coef_diff = 0;
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
model.grid.dealias_method = dealias_method; %de-aliasing method
%model.Smag.dealias_ratio_mask_LS = dealias_ratio_mask_LS;
model.plots = plots_bool;
model.advection.coef_diff = coef_diff;

%% Generating initial buoyancy
[~,model_HR] = fct_buoyancy_init(model,resolution);
[~,model_LR] = fct_buoyancy_init(model,resolution_LR);
clear model

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
if model_HR.advection.HV.bool
    add_subgrid_deter = '_HV';
elseif model_HR.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
    %add_subgrid_deter = [];
end
% if ( model_HR.advection.HV.bool || model_HR.advection.Lap_visco.bool) && ...
%         model_HR.advection.Smag.bool
if ( model_HR.advection.HV.bool | model_HR.advection.Lap_visco.bool) & ...
        model_HR.advection.Smag.bool
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
    %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
    %         fct_num2str(model_HR.advection.Smag.kappamax_on_kappad) ...
    %         '_dealias_ratio_mask_LS_' ...
    %         fct_num2str(model_HR.grid.dealias_ratio_mask_LS)];
end

if isinf(model_HR.sigma.k_c) % Deterministic case
    model_HR.folder.folder_simu = [ 'images/usual_' model_HR.dynamics ...
        add_subgrid_deter '/' model_HR.type_data ];
else % Stochastic case
    model_HR.folder.folder_simu = [ 'images/' model_HR.dynamics ...
        '_MU' add_subgrid_deter '/' model_HR.type_data ];
end
if model_HR.advection.forcing.bool
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_forced_turb' ];
else
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_free_turb' ];
end
model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/' num2str(model_HR.grid.MX(1)) 'x' num2str(model_HR.grid.MX(2)) ];
if ( model_HR.advection.HV.bool | model_HR.advection.Lap_visco.bool) & ...
        model_HR.advection.Smag.bool
    subgrid_details = ['kappamax_on_kappad_' ...
        fct_num2str(model_HR.advection.Smag.kappamax_on_kappad) ...
        '_dealias_ratio_mask_LS_' ...
        fct_num2str(model_HR.advection.Smag.dealias_ratio_mask_LS)];
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '/' subgrid_details ];
end

%model.folder.folder_simu_ref = model.folder.folder_simu;
model_LR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/Low_Pass_fitlered_version_' ...
    num2str(model_LR.grid.MX(1)) 'x' num2str(model_LR.grid.MX(2)) ];

% Create the folders
fct_create_folder_plots(model_LR)

% Colormap
load('BuYlRd.mat');
model_LR.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model_LR.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid

% Spatial grid
model_LR.grid.x = model_LR.grid.dX(1)*(0:model_LR.grid.MX(1)-1);
model_LR.grid.y = model_LR.grid.dX(2)*(0:model_LR.grid.MX(2)-1);

% Grid in Fourier space
model_LR = init_grid_k (model_LR);

% Ensemble size
N_ech=model_LR.advection.N_ech;

%%

t_last_plot = -inf;
%model.advection.step='finite_variation';

t_ini=1;

% model_ref = model;
% name_file = [model.folder.folder_simu_ref '/files/' num2str(0) '.mat'];
name_file_HR = [model_HR.folder.folder_simu '/files/' num2str(0) '.mat'];
load(name_file_HR)
model_HR = model; clear model;
% model.folder.folder_simu_ref = model_ref.folder.folder_simu_ref;

dt=model_HR.advection.dt_adv;
N_t = ceil(model_HR.advection.advection_duration/dt);

%x = model_LR.grid.dX(1)*(0:model.grid.MX(1)-1);
% % if model.mirror
% %     y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
% % else
% %     y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% % end
% My=model.grid.MX(2);
% % if model.mirror
% %     My=My/2;
% % end
% 
% F_save = [];
% F_save2 = [];
% bt1_HR_vect = [];
% bt1_LR_vect = [];
% 
% 
% trigger = false;
% %t_ini=1700000
for t_loop=t_ini:N_t
    %     %% Plot
    % t_loop=1;
    if (t_loop - t_last_plot)*dt >= 3600*24*1
        day = num2str(floor(t_loop*dt/24/3600));
        fprintf([ num2str(t_loop*dt/(24*3600)) ' days of advection \n'])
        
%       model_LR.advection.plot_modes = plot_modes;
%        model_LR.advection.nb_modes = nb_modes;
        t_last_plot=t_loop;
%         id_part=1;
%         
%         width=1.2e3;
%         height=0.5e3;
%         if strcmp(model.type_data,'klein')
%             width=width/2;
%             r_c_ax = 0.5;
%         else
%             r_c_ax =1/1.5;
%             %             r_c_ax =1;
%         end
%         %         X0=0.5e3*[1 1];
%         X0=[0 1];
        
        %% Specific day
        %         warning('the day is specified manually');
        %         day = day_choose;
        %         day = num2str(day);
        
        %% Load
%         model_ref = model;
%         name_file = [model.folder.folder_simu_ref '/files/' num2str(day) '.mat'];
        name_file_HR = [model_HR.folder.folder_simu '/files/' num2str(day) '.mat'];
        load(name_file_HR)
        model_HR = model; clear model;
%         model = model_ref;
        
        %model.odg_b = model.odg_b*3;
        %
        %         fct_plot(model,fft_buoy_part,day)
        
        %% Anti-aliasing filtering and subsampling
        % Filtering
        fft_buoy_part((model_LR.grid.MX(1)/2+1): ...
            (end-1-model_LR.grid.MX(1)/2+1),:,:,:)=0;
        fft_buoy_part(:,(model_LR.grid.MX(2)/2+1): ...
            (end-1-model_LR.grid.MX(2)/2+1),:,:)=0;
        % Subsampling
        buoy_part = real(ifft2(fft_buoy_part));
        clear fft_buoy_part
        buoy_part = buoy_part( ...
            1:model_HR.grid.MX(1)/model_LR.grid.MX(1):end, ...
            1:model_HR.grid.MX(2)/model_LR.grid.MX(2):end,:,:);
        fft_buoy_part_ref = fft2(buoy_part);    
        clear buoy_part;
              
        %%
        if model_LR.advection.Smag.bool
            [coef_diff_aa,coef_diff] = fct_coef_diff(model_LR,fft_buoy_part_ref);
            figure(9);
            subplot(1,2,1);
            imagesc(model_LR.grid.x_ref,model_LR.grid.y_ref,...
                real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
            subplot(1,2,2);
            imagesc(model_LR.grid.x_ref,model_LR.grid.y_ref,...
                coef_diff_aa');axis xy;axis equal;colorbar;
            drawnow
            eval( ['print -depsc ' model_LR.folder.folder_simu ...
                '/dissip_coef/' day '.eps']);
        end
        %%
        % Plots
        spectrum_ref = fct_plot(model_LR,fft_buoy_part_ref,day);
        
%         if model_LR.advection.plot_dissip
%             fct_plot_dissipation(model_LR,fft_buoy_part_ref,sigma_on_sq_dt,day);
%         end
%         
% %         %% Distance to reference
% %         % Spectrum discrepancy
% %         spectrum_ref = ref{1};
% %         error_vs_t(1,1) = dkappa * ...
% %             sum(abs(spectrum(:) - spectrum_ref(:)));
% %         % RMSE
% %         fft_buoy_ref = ref{2};
% %         error_vs_t(1,2) = 1/prod(model.grid.MX)^2 * ...
% %             sum(abs(fft_buoy_part_ref(:) - fft_buoy_ref(:)).^2);
% %         % Square root
% %          error_vs_t = sqrt( error_vs_t);
         
        %% Save files
        save( [model_LR.folder.folder_simu '/files/' day '.mat'], ...
            'model_LR','t','fft_buoy_part_ref','spectrum_ref');
        
    end
end
end