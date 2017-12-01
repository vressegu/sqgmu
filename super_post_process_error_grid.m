
%%%%%%%%%%%%%%%%%%%%
%%% Super main
%%%%%%%%%%%%%%%%%%%%
init;

%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
% dynamics = '2D';

nb_days =30;

% Deterministic or random model
stochastic_simulation = false;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

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
resolution = 128;
%resolution = 1024;
%resolution = 2048;

resolution_HR = 1024;

% The number of grid point is resolution^2
% It has to be an even integer

% Forcing

% Forcing or not
forcing = false;
% If yes, there is a forcing
% F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y)
% % If yes, there is an additionnal velocity V = (0 Vy)
% % with Vy = ampli_forcing * odg_b *  sin( 2 freq_f pi y/L_y)

% Viscosity
Lap_visco.bool = true;

% % Smagorinsky-like viscosity
% Smag.bool = false;
% % HV.bool = false;

% Hyper-viscosity
HV.bool = false;

% Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
Smag.bool = true;

if Lap_visco.bool
    % Ratio between the Shanon resolution and filtering frequency used to
    % filter the heterogenous diffusion coefficient
    v_dealias_ratio_mask_LS = 1./ [1 2 4 8 16 64]';
    v_dealias_ratio_mask_LS = sort(v_dealias_ratio_mask_LS);
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
    % and the targeted diffusion scale
    v_kappamax_on_kappad = 0.4:0.2:1.2 ;
    v_kappamax_on_kappad = sort(v_kappamax_on_kappad);
    
    for p=1:length(v_dealias_ratio_mask_LS)
        for q=1:length(v_kappamax_on_kappad)
            Smag(p,q).bool = true;
            Smag(p,q).dealias_ratio_mask_LS = v_dealias_ratio_mask_LS(p);
            Smag(p,q).kappamax_on_kappad = v_kappamax_on_kappad(q);
            Smag(p,q).weight_cst_dissip = 0;
        end
    end
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


%% Main
s_Smag = size(Smag);
Smag = Smag(:);
%Smag = Smag(1:2);
ll=length(Smag);
error_vs_t = nan([nb_days,2,ll]);
for j=1:ll
%parfor j=1:ll
    error_vs_t(:,:,j) = post_process_error_grid(...
        stochastic_simulation,type_data,resolution,resolution_HR,...
        forcing,sigma(j),Lap_visco,HV,Smag(j));
end
error_vs_t = reshape(error_vs_t,[nb_days 2 s_Smag]);

%% Folder to save plots and files


% Physical parameters
model = fct_physical_param(dynamics);

% Gather parameters in the structure model
% model.sigma.slope_sigma = slope_sigma;
% model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
% if strcmp(type_data,'Spectrum')
%     model.slope_b_ini = slope_b_ini;
% end
model.dynamics=dynamics;
model.type_data=type_data;
% model.advection.N_ech=N_ech;
% model.sigma.k_c = k_c;
% model.advection.advection_duration=advection_duration;
% model.advection.plot_dissip = plot_dissip;
% model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
% model.advection.forcing.ampli_forcing = ampli_forcing;
% model.advection.forcing.freq_f = freq_f;
% model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag(1);
% model.advection.use_save = use_save;
% model.advection.day_save = day_save;
% model.grid.dealias_method = dealias_method; %de-aliasing method
% %model.Smag.dealias_ratio_mask_LS = dealias_ratio_mask_LS;
% model.plots = plots_bool;
% model.advection.coef_diff = coef_diff;


if model.advection.HV.bool
    add_subgrid_deter = '_HV';
elseif model.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
    %add_subgrid_deter = [];
end
% if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
%         model.advection.Smag.bool
if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
    %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
    %         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
    %         '_dealias_ratio_mask_LS_' ...
    %         fct_num2str(model.grid.dealias_ratio_mask_LS)];
end

if ~stochastic_simulation % Deterministic case
    model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
else % Stochastic case
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_MU' add_subgrid_deter '/' model.type_data ];
end
if model.advection.forcing.bool
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_forced_turb' ];
else
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_free_turb' ];
end
model.folder.folder_simu = [ model.folder.folder_simu ...
    '/' num2str(resolution) 'x' num2str(resolution) ];
% if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%         model.advection.Smag.bool
%     subgrid_details = ['kappamax_on_kappad_' ...
%         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
%         '_dealias_ratio_mask_LS_' ...
%         fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
%     model.folder.folder_simu = [ model.folder.folder_simu ...
%         '/' subgrid_details ];
% end
% % Create the folders
% fct_create_folder_plots(model)

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

%%

% if HV.bool
%     add_subgrid_deter = '_HV';
% elseif model.advection.Lap_visco.bool
%     add_subgrid_deter = '_Lap_visco';
% else
%     add_subgrid_deter = '_no_deter_subgrid';
%     %add_subgrid_deter = [];
% end
% % if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
% %         model.advection.Smag.bool
% if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%         model.advection.Smag.bool
%     add_subgrid_deter = [add_subgrid_deter '_Smag'];
%     %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
%     %         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
%     %         '_dealias_ratio_mask_LS_' ...
%     %         fct_num2str(model.grid.dealias_ratio_mask_LS)];
% end
% 
% if ~stochastic_simulation % Deterministic case
%    folder_simu = [ 'images/usual_' dynamics ...
%         add_subgrid_deter '/' type_data ];
% else % Stochastic case
%     folder_simu = [ 'images/' dynamics ...
%         '_MU' add_subgrid_deter '/' type_data ];
% end
% if forcing.bool
%     folder_simu = [ folder_simu ...
%         '_forced_turb' ];
% else
%     folder_simu = [folder_simu ...
%         '_free_turb' ];
% end
% folder_simu = [ folder_simu ...
%     '/' num2str(resolution) 'x' num2str(resolution) ];
% % if ( HV.bool | model.advection.Lap_visco.bool) & ...
% %         model.advection.Smag.bool
% %     subgrid_details = ['kappamax_on_kappad_' ...
% %         fct_num2str(Smag.kappamax_on_kappad) ...
% %         '_dealias_ratio_mask_LS_' ...
% %         fct_num2str(Smag.dealias_ratio_mask_LS)];
% %     folder_simu = [ folder_simu ...
% %         '/' subgrid_details ];
% % end
% % Create the folders
% model.folder.folder_simu = folder_simu;
% fct_create_folder_plots(model)


%% Save
mean_error = squeeze(mean(error_vs_t,1));
mkdir([model.folder.folder_simu '/error_files']);
save( [model.folder.folder_simu '/error_files/error.mat'], ...
    'model','Smag','v_dealias_ratio_mask_LS','v_kappamax_on_kappad',...
    'error_vs_t','mean_error');

%% Plots
close all
taille_police = 12;

figure(1)
imagesc(v_dealias_ratio_mask_LS,v_kappamax_on_kappad,...
    squeeze(mean_error(1,:,:))');
axis xy;axis equal; colorbar;
ylabel('$(\pi/\Delta x)/k_d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('Filtering diff. coef.',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('time-mean Spectrum error (m.s$^-2$)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')

figure(2)
imagesc(v_dealias_ratio_mask_LS,v_kappamax_on_kappad,...
    squeeze(mean_error(2,:,:))');
axis xy;axis equal; colorbar;
ylabel('$(\pi/\Delta x)/k_d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('Filtering diff. coef.',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('time-mean RMSE (m.s$^-2$)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')



%% Plots at 7 days
error_7days = squeeze(error_vs_t(8,:,:,:));

figure(3)
imagesc(v_dealias_ratio_mask_LS,v_kappamax_on_kappad,...
    squeeze(error_7days(1,:,:))');
axis xy;axis equal; colorbar;
ylabel('$(\pi/\Delta x)/k_d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('Filtering diff. coef.',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('Spectrum error (m.s$^-2$) at $7$ days',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')

figure(4)
imagesc(v_dealias_ratio_mask_LS,v_kappamax_on_kappad,...
    squeeze(error_7days(2,:,:))');
axis xy;axis equal; colorbar;
ylabel('$(\pi/\Delta x)/k_d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('Filtering diff. coef.',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('RMSE (m.s$^-2$) at $7$ days',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')



