function [fft_b, model] = fct_fft_advection_sto_mat(model,  fft_b)
% Advection of buoyancy using SQG or SQG_MU model
%

tic
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
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_V2
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_energy_flux
    add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
    if model.sigma.hetero_energy_flux_v2
        add_subgrid_deter = [add_subgrid_deter '_v2'];
    end
    if model.sigma.hetero_energy_flux_averaging_after
        add_subgrid_deter = [add_subgrid_deter '_1_3rd_before_norm'];
    end
    if isfield(model.sigma,'kappa_VLS_on_kappa_LS')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappa_LS_on_kappa_VLS_' ...
            num2str(1/model.sigma.kappa_VLS_on_kappa_LS)];
    end
    if isfield(model.sigma,'kappaLSforEspi_on_kappamin')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappamin_on_kappaLSforEspi__' ...
            num2str(1/model.sigma.kappaLSforEspi_on_kappamin)];
    end
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_Smag
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_Smag'];
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
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
    if ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF'))
        subgrid_details = [ subgrid_details ...
            '_kappamin_on_kappamax_' ....
            fct_num2str(model.sigma.kappamin_on_kappamax) ];
        if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
            subgrid_details = [ subgrid_details ...
                '_on_kc_' ....
                fct_num2str(1/model.sigma.k_c) ];
        elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                if model.sigma.estim_k_LS
                    subgrid_details = [ subgrid_details ...
                     '_estim_k_LS'];
                end
                if model.sigma.time_smooth.bool
                    subgrid_details = [ subgrid_details ...
                        '_time_smooth_'... 
                        num2str(24*3600/model.sigma.time_smooth.tau)];
                end
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
    if model.advection.N_ech > 1
        subgrid_details = [ subgrid_details ...
            '_N_ech_' ....
            fct_num2str(model.advection.N_ech) ];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end

% Create the folders
fct_create_folder_plots(model)

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

%% Initialisation of the spatial fields

% Remove aliasing
fft_b(model.grid.k.ZM(1),:,:,:)=0;
fft_b(:,model.grid.k.ZM(2),:,:)=0;

% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_b);
w=real(ifft2(fft_w));

% figure;imagesc(real(ifft2(fft_b))');axis xy;axis equal;colorbar;
% figure;imagesc(w(:,:,1)');axis xy;axis equal;colorbar;

% Create several identical realizations of the intial buoyancy
fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);

%% Choice of the variance tensor a
if model.sigma.sto & ...
        ( model.sigma.Smag.bool | model.sigma.assoc_diff )
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        [~, ~, tr_a ,...
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ...
            ]= fct_sigma_spectrum_abs_diff(...
            model,fft_w,true,num2str(0));
        %         a0 = tr_a/2;
        %         model.sigma.a0 = max(a0);
        a0_LS_temp = tr_a/2;
    else
        % Variance tensor of the simulated small-scle velocity
        a0_LS_temp = 1/2 * model.sigma.fct_tr_a(model, ...
            model.sigma.kappamin_on_kappamax ...
            * model.sigma.kappamax_on_kappaShanon ...
            * (pi/sqrt(prod(model.grid.dX))), ...
            model.sigma.kappamax_on_kappaShanon ...
            * (pi/sqrt(prod(model.grid.dX))));
        %     a0_LS_temp = 1/2 * fct_norm_tr_a_theo(model, ...
        %         model.sigma.kappamin_on_kappamax ...
        %         * model.sigma.kappamax_on_kappaShanon ...
        %         * (pi/sqrt(prod(model.grid.dX))), ...
        %         model.sigma.kappamax_on_kappaShanon ...
        %         * (pi/sqrt(prod(model.grid.dX))), ...
        %         ( 3 - model.sigma.slope_sigma )/2 );
        % Variance tensor of the unsimulated small-scle velocity
        %     a0_SS_temp = 1/2 * fct_norm_tr_a_theo(model, ...
    end
    a0_SS_temp = 1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(model, ...
        model.sigma.kappaMinUnresolved_on_kappaShanon ...
        * (pi/sqrt(prod(model.grid.dX))), ...
        model.sigma.kappaMaxUnresolved_on_kappaShanon ...
        * (pi/sqrt(prod(model.grid.dX))));
    % Multiplicative foctor for the dissipation coefficient
    if model.sigma.assoc_diff
        % Variance tensor
        if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            model.sigma.a0 = a0_LS_temp + a0_SS_temp ;
        else
            model.sigma.a0 = ( 1 + a0_SS_temp / a0_LS_temp )...
                * 2 * model.physical_constant.f0 / model.sigma.k_c^2;
        end
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
    elseif model.sigma.Smag.bool
        if model.sigma.Smag.epsi_without_noise
            coef_diff_temp = 1;
        else
            if a0_SS_temp > eps
                coef_diff_temp = 1 + a0_LS_temp / a0_SS_temp ;
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                coef_diff_temp = 0;
            else
                error('Unknow case');
            end
        end
        % Heterogeneous dissipation coefficient
        coef_diff_aa_temp = coef_diff_temp * ...
            fct_coef_diff(model,fft_b);
        % Coefficient coef_Smag to target a specific diffusive scale
        model.sigma.Smag.coef_Smag = ...
            ( model.sigma.Smag.kappamax_on_kappad ...
            * sqrt(prod(model.grid.dX))/pi ) ^ 2;
        coef_diff_aa_temp = model.sigma.Smag.coef_Smag * ...
            coef_diff_aa_temp ;
        % Maximum of the total variance tensor
        model.sigma.a0 = 2 * max(coef_diff_aa_temp(:));
        clear a0_LS_temp a0_SS_temp coef_diff_temp coef_diff_aa_temp
    else
        error('Unknown case');
    end
    
elseif model.sigma.sto & ...
        strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
    [~, ~, tr_a ,....
        model.sigma.slope_sigma,...
        model.sigma.offset_spectrum_a_sigma, ...
        model.sigma.km_LS ...
        ] = fct_sigma_spectrum_abs_diff(model,fft_w,true,num2str(0));
    % sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) * sigma; clear sigma
    a0 = tr_a/2;
    model.sigma.a0 = max(a0);
    %     model.sigma.k_c = ...
    %         sqrt( 2 * model.physical_constant.f0 / model.sigma.a0 );
    %     % Diffusion coefficient
    %     model.advection.coef_diff = 1/2 * model.sigma.a0;
elseif model.sigma.sto & ...
        ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF'))
    % sigma = nan;
    
    current_folder = pwd;
    cd(model.folder.folder_simu);
    cd ..
    model.folder.folder_EOF = [ pwd '/folder_EOF'];
    cd(current_folder);
    
    % Load precomputed EOFs and correspond variance tensor
    load([ model.folder.folder_EOF '/EOF.mat'],'EOF');
    
    %%
    
%     load([ model.folder.folder_EOF '/var_tensor_from_EOF.mat'],...
%         'a_xx','a0');
%     
% %     a_xx = sum( permute(EOF,[1 2 3 5 4]) .* permute(EOF,[1 2 5 3 4]) ,5);
% %     a02 =  1/2 * ( a_xx(:,:,1,1) + a_xx(:,:,2,2) );
% %     a02 =  mean(a02(:));
    
    EOF = permute(EOF,[1 2 3 5 4]);
    if isinf(model.sigma.nb_EOF)
        model.sigma.nb_EOF = size(EOF,5);
    else
        EOF = EOF(:,:,:,:,1:model.sigma.nb_EOF);
    end    
    %model.sigma.nb_EOF = size(EOF,5);
    sigma = EOF; clear EOF;
    
    a_xx = sum( sigma .* permute(sigma,[1 2 4 3 5]) ,5);
    a0 =  1/2 * ( a_xx(:,:,1,1) + a_xx(:,:,2,2) );
    a0 =  mean(a0(:));
    
    % Filtering the variance tensor at large scales
    sigma_loc = permute(a_xx,[3 4 1 2]);
    sigma_loc = multi_chol(sigma_loc);
    sigma_loc = permute(sigma_loc,[3 4 1 2]);
    sigma_loc_aa = fft2(sigma_loc);
    sigma_loc_aa = real(ifft2( bsxfun(@times, sigma_loc_aa,...
        model.grid.k_aa_LS.mask) ));
    sigma_loc_aa = permute(sigma_loc_aa,[3 4 1 2]);
    a_xx_aa = multiprod(sigma_loc_aa,multitrans(sigma_loc_aa));
    a_xx_aa = permute(a_xx_aa,[3 4 1 2]);
    
    %     % Set negative values to zero
    %     a_xx = fct_smooth_pos_part(a_xx);
    %     a_xx_aa = fct_smooth_pos_part(a_xx_aa);
    
    % Plots
    figure(12);
    for d1=1:2
        for d2=1:2
            subplot(2,2,1+d1-1+2*(d2-1))
            imagesc(model.grid.x,model.grid.y,a_xx(:,:,d1,d2)');
            axis equal; axis xy;colorbar;
            title(['$a_{' num2str(d1) num2str(d2) '}$'],...
                'FontUnits','points',...
                'FontWeight','normal',...
                'interpreter','latex',...
                'FontSize',12,...
                'FontName','Times')
        end
    end
    drawnow
    eval( ['print -depsc ' model.folder.folder_simu ...
        '/Variance_tensor.eps']);
    figure(13);
    for d1=1:2
        for d2=1:2
            subplot(2,2,1+d1-1+2*(d2-1))
            imagesc(model.grid.x,model.grid.y,a_xx_aa(:,:,d1,d2)');
            axis equal; axis xy;colorbar;
            title(['Smooth $a_{' num2str(d1) num2str(d2) '}$'],...
                'FontUnits','points',...
                'FontWeight','normal',...
                'interpreter','latex',...
                'FontSize',12,...
                'FontName','Times')
        end
    end
    drawnow
    eval( ['print -depsc ' model.folder.folder_simu ...
        '/Smooth_Variance_tensor.eps']);
    
    % Variance tensor
    [a0 max(a_xx_aa(:)) ]
    a0 = max( [a0 max(a_xx_aa(:)) ]);
    model.sigma.a0 = a0; clear a0;
    % Diffusion coefficient
    a_xx_aa = permute( a_xx_aa, [ 1 2 4 5 3]);
    model.advection.coef_diff = 1/2 * a_xx_aa; clear a_xx;
else
    % Variance tensor
    model.sigma.a0 = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
    % Diffusion coefficient
    model.advection.coef_diff = 1/2 * model.sigma.a0;
end


if model.sigma.sto & model.sigma.hetero_energy_flux
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        model.sigma.km_LS = repmat(model.sigma.km_LS,[1 1 1 N_ech]);
    end
    coef_modulation = fct_epsilon_k_onLine(model,fft_b);
elseif model.sigma.sto & model.sigma.hetero_modulation_Smag
    % Heterogeneous dissipation coefficient
    coef_modulation = fct_coef_diff(model,fft_b);
    m_coef_modulation = mean(mean(coef_modulation,1),2);
    coef_modulation = bsxfun( @times, ...
        1./m_coef_modulation, coef_modulation);
    clear m_coef_modulation
    % coef_modulation = fct_coef_diff(model,fft_b);
elseif model.sigma.sto & ...
        (model.sigma.hetero_modulation | model.sigma.hetero_modulation_V2)
    coef_modulation = fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
else
    coef_modulation = 1;
end
model.sigma.a0 = model.sigma.a0 * max(coef_modulation(:));
% warning(['The CFL may changed here ...
% in the case of heterogreneous variance tensor'])

%% Forcing
% if isfield(model.advection, 'forcing') && model.advection.forcing.bool
%     %     Ly = model.grid.MX(2) * model.grid.dX(2);
%     %     [~,Y]=ndgrid(model.grid.x,model.grid.y);
%     %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
%     %         sin( 2 * pi * model.advection.forcing.freq_f * Y / Ly);
%     %     clear Ly X Y
%
%
%     model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
%     model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
%     [X,Y]=ndgrid(model.grid.x,model.grid.y);
%
%     switch model.dynamics
%         case 'SQG'
%             U_caract = ...
%                 model.odg_b / model.physical_constant.buoyancy_freq_N;
%         case '2D'
%             U_caract = ...
%                 model.odg_b * model.advection.forcing.Ly;
%             U_caract = U_caract /20;
%             U_caract = U_caract /4;
%         otherwise
%             error('Unknown type of dynamics');
%     end
%     %     U_caract = U_caract /10;
%     U_caract = U_caract /3;
%
%     %     model.advection.U_caract = sqrt(mean(w(:).^2));
%     model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
%
%     %     model.advection.forcing.on_T = model.advection.forcing.on_T /10;
%
%     model.advection.forcing.ampli_forcing = ...
%         model.advection.forcing.ampli_forcing * model.odg_b ;
%     %     model.advection.forcing.ampli_forcing = ...
%     %         model.advection.forcing.ampli_forcing * model.odg_b * ...
%     %         model.advection.forcing.on_T;
%
%     %     ampli_scale = sqrt( ...
%     %         ( sum(model.advection.forcing.freq_f.^2)/2 ) ...
%     %             ^(model.sigma.slope_sigma) ...
%     %             );
%     ampli_scale = 1;
%
%     switch model.advection.forcing.forcing_type
%         case 'Kolmogorov'
%             F = ampli_scale * ...
%                 model.advection.forcing.ampli_forcing * ...
%                 cos( 2 * pi / model.advection.forcing.Lx ...
%                 * model.advection.forcing.freq_f(1) * X ...
%                 + 2 * pi / model.advection.forcing.Ly ...
%                 * model.advection.forcing.freq_f(2) * Y );
%             F = F / model.advection.forcing.on_T;
%         case 'Spring'
%             F = ampli_scale * ...
%                 model.advection.forcing.ampli_forcing ...
%                 * sin( 2 * pi / model.advection.forcing.Lx ...
%                 * model.advection.forcing.freq_f(1) * X ) ...
%                 .* sin( 2 * pi / model.advection.forcing.Ly ...
%                 * model.advection.forcing.freq_f(2) * Y );
%             %     F = model.advection.forcing.ampli_forcing * ...
%             %         sin( 2 * pi * model.advection.forcing.freq_f * Y / ...
%             %         model.advection.forcing.Ly);
%     end
%
%
%
%     model.advection.forcing.F = fft2(F);
%
%     if strcmp(model.type_data, 'Zero')
%         fft_w = SQG_large_UQ(model,  ...
%             model.odg_b / model.advection.forcing.ampli_forcing ...
%             * model.advection.forcing.F);
%         %   w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
%         %         w(:,:,2) = 0;
%         %         fft_w = fft2(w);
%         if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
%             fft_w = fft_w * model.advection.forcing.on_T;
%         end
%         w = real(ifft2(fft_w));
%     end
%
%     %
%     %     figure;imagesc(model.grid.x,model.grid.y,model.advection.forcing.F');
%     %     axis xy; axis equal;
%     %
%
%     %     model.advection.forcing.F = fft2(F);
%     clear Lx Ly X Y on_T U_caract ampli_scale
% end
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %     Ly = model.grid.MX(2) * model.grid.dX(2);
    %     [~,Y_forcing]=ndgrid(model.grid.x,model.grid.y);
    %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         sin( 2 * pi * model.advection.forcing.freq_f * Y_forcing / Ly);
    %     clear Ly X_forcing Y_forcing
    
    
    model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
    model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
    [X_forcing,Y_forcing]=ndgrid(model.grid.x,model.grid.y);
    
    
    
    switch model.dynamics
        case 'SQG'
            U_caract = ...
                model.odg_b / model.physical_constant.buoyancy_freq_N;
            U_caract = U_caract /5;
        case '2D'
            U_caract = ...
                model.odg_b * model.advection.forcing.Ly;
            U_caract = U_caract /20;
            U_caract = U_caract /4;
        otherwise
            error('Unknown type of dynamics');
    end
    %     U_caract = U_caract /10;
    U_caract = U_caract /3;
    
    %     model.advection.U_caract = sqrt(mean(w(:).^2));
    model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
    
    %     model.advection.forcing.on_T = model.advection.forcing.on_T /10;
    
    model.advection.forcing.ampli_forcing = ...
        model.advection.forcing.ampli_forcing * model.odg_b ;
    %     model.advection.forcing.ampli_forcing = ...
    %         model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         model.advection.forcing.on_T;
    
    %     ampli_scale = sqrt( ...
    %         ( sum(model.advection.forcing.freq_f.^2)/2 ) ...
    %             ^(model.sigma.slope_sigma) ...
    %             );
    ampli_scale = 1;
    
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing * ...
                cos( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ...
                + 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
            F = F / model.advection.forcing.on_T;
        case {'Spring','Hetero_Spring'}
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing ...
                * sin( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ) ...
                .* sin( 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
            %     F = model.advection.forcing.ampli_forcing * ...
            %         sin( 2 * pi * model.advection.forcing.freq_f * Y_forcing / ...
            %         model.advection.forcing.Ly);
    end
    
    
    if strcmp(model.advection.forcing.forcing_type, 'Hetero_Spring')
        model.advection.forcing.F = F;
    else
        model.advection.forcing.F = fft2(F);
    end
    
    if strcmp(model.type_data, 'Zero')
        fft_w = SQG_large_UQ(model,  ...
            model.odg_b / model.advection.forcing.ampli_forcing ...
            * model.advection.forcing.F);
        %         w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
        %         w(:,:,2) = 0;
        %         fft_w = fft2(w);
        if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
            fft_w = fft_w * model.advection.forcing.on_T;
        end
        w = real(ifft2(fft_w));
    end
    
    %
    %     figure;imagesc(model.grid.x,model.grid.y,model.advection.forcing.F');
    %     axis xy; axis equal;
    %
    
    %     model.advection.forcing.F = fft2(F);
    clear Lx Ly X_forcing Y_forcing on_T U_caract ampli_scale
end

%% Hyperviscosity
if model.advection.Lap_visco.bool | model.advection.HV.bool
    % Root mean square Lyapunov
    [dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
    [dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
    s=(dxUx+dyUy)/2;
    d=(dyUx+dxUy)/2;
    lambda = sqrt(s.^2+d.^2);
    model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
    clear s d dxUx dxUy dyUx dyUy lambda
    
    % Order of the hyperviscosity
    if model.advection.HV.bool
        % model.advection.HV.order=8;
    else
        model.advection.HV.order=2;
    end
    
    % Hyperviscosity coefficient
    model.advection.HV.val= ...
        40 * model.advection.lambda_RMS * ...
        (mean(model.grid.dX)/pi)^model.advection.HV.order;
    
    model.advection.HV.val
    % model.advection.HV.val*(model.grid.MX(1)^model.advection.HV.order)
    model.advection.HV.val/(mean(model.grid.dX)^model.advection.HV.order)
    model.advection.lambda_RMS
    
    % if ~model.advection.HV.bool % DNS
    if model.advection.Lap_visco.bool &&  strcmp(model.dynamics,'2D')
        if isfield(model.advection, 'forcing') && model.advection.forcing.bool
            model.carct.L_caract = 1/ ...
                sqrt(mean( ( [ model.advection.forcing.Lx ...
                model.advection.forcing.Ly ] ...
                .\ model.advection.forcing.freq_f ).^2));
            warning('Re should be computed with the amplitude of forcing!');
        else
            fft_grad_b = fct_grad(model,fft_b);
            model.carct.L_caract = sqrt( ...
                sum(abs(fft_b(:)).^2) / sum(abs(fft_grad_b(:)).^2) );
        end
        model.carct.U_caract = model.odg_b * model.carct.L_caract;
        Re = model.carct.U_caract * model.carct.L_caract ...
            / model.advection.HV.val
        model.carct.Re = Re;
        model.carct.l_Kolmogorov = Re^(-1/2) * model.carct.L_caract;
        if model.carct.l_Kolmogorov < min(model.grid.dX)
            error('The simulation is under resolved');
        end
    end
    
    if model.advection.Smag.bool
        model.advection.Smag.coef_Smag = ...
            ( model.advection.Smag.kappamax_on_kappad ...
            * sqrt(prod(model.grid.dX))/pi ) ...
            ^ (5*model.advection.HV.order/4 -1/2);
        
        
        % Heterogeneous HV or diffusivity/viscosity coefficient
        if model.advection.Lap_visco.bool
            % Heterogeneous dissipation coefficient
            coef_diff_aa = fct_coef_diff(model, fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            coef_diff_aa = model.advection.Smag.coef_Smag * coef_diff_aa ;
            % Maximum dissipation coefficient
            model.advection.HV.maxVal = max(coef_diff_aa(:));
        elseif model.advection.HV.bool
            % Heterogeneous HV coefficient
            coef_HV_aa = fct_coef_HV(model, fft_b);
            
            % Maximum HV coefficient
            model.advection.HV.maxVal = max(coef_HV_aa(:));
        else
            error('Unknown deterministic subgrid tensor');
        end
    else
        model.advection.HV.maxVal = model.advection.HV.val;
    end
else
    model.advection.HV.val= 0;
    model.advection.HV.order = 2;
    model.advection.HV.maxVal = 0;
end

%% Choice of time step : CFL
model.advection.dt_adv = fct_CFL(model,w);

%% Fourier transform of the kernel \tilde sigma
if model.sigma.sto
    %if model.sigma.SelfSim_from_LS.bool
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        %[sigma, ~, tr_a ] ...
        [sigma, ~, tr_a ,....
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ]...
            = fct_sigma_spectrum_abs_diff(model,fft_w,true,num2str(0));
        % sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) * sigma; clear sigma
        sigma = repmat(sigma,[1 1 1 N_ech]);
        model.sigma.km_LS = repmat(model.sigma.km_LS,[1 N_ech]);
        a0 = tr_a/2;
        model.sigma.a0 = a0;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        missed_var_small_scale_spectrum = 2*a0;
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
    elseif ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') )
        % Fourier transform of the kernel \tilde sigma up to a multiplicative
        % constant
        %     [sigma_on_sq_dt, ~, missed_var_small_scale_spectrum ] ...
        [sigma, ~, missed_var_small_scale_spectrum ] ...
            = fct_sigma_spectrum(model,fft_w);
        %   warning('Above line modified for debug !!!!!')
        % sigma_on_sq_dt will be used to simulate sigma d B_t
        % missed_var_small_scale_spectrum will be used to set the mulitplicative
        % constant
    end
    
    if model.sigma.Smag.bool | model.sigma.assoc_diff
        % Variance tensor of the simulated small-scle velocity
        model.sigma.a0_LS = 1/2 * missed_var_small_scale_spectrum;
        % Variance tensor of the unsimulated small-scle velocity
        %         model.sigma.a0_SS = 1/2 * fct_norm_tr_a_theo(...
        model.sigma.a0_SS = 1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
            model, ...
            model.sigma.kappaMinUnresolved_on_kappaShanon ...
            *(pi/sqrt(prod(model.grid.dX))), ...
            model.sigma.kappaMaxUnresolved_on_kappaShanon ...
            *(pi/sqrt(prod(model.grid.dX))));
        model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        
        % Muliplicative constant of the kernel \tilde sigma
        % sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) * sigma; clear sigma
        
        if model.sigma.assoc_diff
            if ~ strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % model.sigma.a0_LS
                if model.sigma.a0_LS < eps
                    error(['This should not happened and' ...
                        'may casue divison by 0']);
                end
                % Variance tensor
                a0_LS = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
                coef_temp = a0_LS / model.sigma.a0_LS; clear a0_LS;
                
                model.sigma.a0_LS = coef_temp * model.sigma.a0_LS;
                model.sigma.a0_SS = coef_temp * model.sigma.a0_SS;
                model.sigma.a0 = coef_temp * model.sigma.a0;
                model.sigma.a0_on_dt = coef_temp * model.sigma.a0_on_dt;
                
                sigma = sqrt(coef_temp) * sigma;
                % sigma_on_sq_dt = sqrt(coef_temp) * sigma_on_sq_dt;
            end
            % Diffusion coefficient
            model.advection.coef_diff = 1/2 * model.sigma.a0;
            
        elseif model.sigma.Smag.bool
            if model.sigma.a0_SS > eps
                if model.sigma.Smag.epsi_without_noise
                    sigma = ...
                        sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                        * sigma;
                    %                     sigma_on_sq_dt = ...
                    %                         sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                    %                         * sigma_on_sq_dt;
                    model.advection.coef_diff = 1;
                else
                    sigma = sqrt(2/model.sigma.a0_SS) ...
                        * sigma;
                    %                     sigma_on_sq_dt = sqrt(2/model.sigma.a0_SS) ...
                    %                         * sigma_on_sq_dt;
                    model.advection.coef_diff = 1 + ...
                        model.sigma.a0_LS / model.sigma.a0_SS ;
                end
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                sigma = zeros(size(sigma));
                % sigma_on_sq_dt = zeros(size(sigma_on_sq_dt));
                model.advection.coef_diff = 0;
            else
                error('Unknow case');
            end
            
            %             % Diffusion coefficient
            %             if model.sigma.Smag.epsi_without_noise
            %                 model.advection.coef_diff = 1;
            %             elseif model.sigma.a0_SS > eps
            %                 model.advection.coef_diff = 1 + ...
            %                     model.sigma.a0_LS / model.sigma.a0_SS ;
            % elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            %                 model.advection.coef_diff = 0;
            %             else
            %                 error('Unknow case');
            %            end
        else
            error('Unknown case');
        end
        %        warning('The product (nu_Smag * sigma dBt)
        %          may need to be regularized')
        %     elseif model.sigma.assoc_diff
        %         fsavsd
        %
        %         model.sigma.a0 = a0;
        %         model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        %
        %         % Diffusion coefficient
        %         model.advection.coef_diff = 1/2 * model.sigma.a0;
    elseif ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') )
        % Muliplicative constant of the kernel \tilde sigma
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        sigma = ...
            sqrt(2*model.sigma.a0/missed_var_small_scale_spectrum) ...
            * sigma;
        %         sigma_on_sq_dt = ...
        %             sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
        %             * sigma; clear sigma
        %         %             * sigma_on_sq_dt;
        % the factor d=2 is for the dimension d of the space R^d
        % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
        
    end
    clear missed_var_small_scale_spectrum
    % warning('This formula may changed if k inf > la resolution');
else
    sigma = 0;
    % sigma_on_sq_dt = 0;
end

%% Possibly reset velocity to zero
% if strcmp(model.type_data, 'Zero')
%     w = zeros(size(w));
% end
% clear w fft_w

%% Loop on time

% Used for the first plot
tt_last = -inf;

% Number of time step
% N_t = ceil(model.advection.advection_duration/model.advection.dt_adv);

%% Print informations

if model.plots
    % Printing some information
    fprintf(['The initial condition is ' model.type_data ' \n'])
    fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
    fprintf(['Time step : ' num2str(model.advection.dt_adv) ' seconds \n']);
    fprintf(['Time of advection : ' num2str(...
        model.advection.advection_duration/3600/24) ' days \n']);
    fprintf(['Ensemble size : ' num2str(N_ech) ' realizations \n']);
    fprintf(['Resolution : ' num2str(model.grid.MX(1)) ' x ' ...
        num2str(model.grid.MX(2)) ' \n']);
    if model.advection.HV.bool
        str_subgridtensor = 'Hyper-viscosity';
    elseif model.advection.Lap_visco.bool
        str_subgridtensor = 'Laplacian diffusion';
    else
        str_subgridtensor = 'None';
    end
    if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
            model.advection.Smag.bool
        str_subgridtensor = ['Heterogneous ' str_subgridtensor];
    end
    fprintf(['Deterministic subgrid tensor : ' str_subgridtensor ' \n']);
    fprintf(['Model type : ' add_subgrid_deter ' \n']);
    if model.sigma.sto | model.advection.Smag.bool
        fprintf(['Details : ' subgrid_details ' \n']);
    end
    if model.sigma.sto
        fprintf(['type spectrum sigma :' model.sigma.type_spectrum ' \n']);
    end
end

if model.sigma.sto & ...
        strcmp(model.sigma.type_spectrum,'SelfSim_from_LS') & ...
        model.sigma.time_smooth.bool
    sigma_s = 0;
    % sigma_s = sigma;
end

%% Use a saved files of a former simulation ?
if model.advection.use_save
    warning(['The run begin from an older file instead of from the' ...
        'initial condition']);
    day = num2str(model.advection.day_save);
    name_file = [model.folder.folder_simu '/files/' day '.mat']; clear day
    clear fft_b
    model_ref = model;
    
    load(name_file)
    %      model.advection.advection_duration =  ...
    %          model_ref.advection.advection_duration;
    %      model = catstruct(model_ref,model);
    model = model_ref;
    
    if (exist('fft_b','var')==1)
    elseif (exist('fft_buoy_part','var')==1)
        fft_b = fft_buoy_part; clear fft_buoy_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_tracer_part','var')==1)
        fft_b = fft_tracer_part; clear fft_tracer_part
    elseif (exist('fft_buoy_part_ref','var')==1)
        fft_b = fft_buoy_part_ref; clear fft_buoy_part_ref
    else
        error('Cannot find buoyancy field')
    end
    if model.sigma.sto
        if size(fft_b,4) < model.advection.N_ech
            if size(fft_b,4) == 1
                fft_b = repmat ( fft_b, [ 1 1 1 model.advection.N_ech ]);
                clear w fft_w
            else
                error('The number of realisation of the saved file is too low');
            end
        end
        if size(fft_b,4) > model.advection.N_ech
            warning(['The number of realisation of the saved file is too high.' ...
                ' Some realisations are hence removed.']);
            fft_b(:,:,:,(model.advection.N_ech+1):end) = [];
        end
        if ( strcmp(model.sigma.type_spectrum,'EOF') || ...
                strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ...
                && ( model.sigma.nb_EOF > model.advection.N_ech )
            warning(['The number of EOF is larger than the ensemble size.' ...
                ' Some EOFs are hence removed.']);
            model.sigma.nb_EOF = model.advection.N_ech;
            sigma(:,:,:,:,(model.advection.N_ech+1):end) = [];
        end
    end
    
    % Version of matlab
    vers = version;
    year = str2double(vers(end-5:end-2));
    subvers = vers(end-1);
    model.folder.colormap_freeze = ...
        (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
    
    if ~isfield(model,'plots')
        model.plots = true;
    end
    
    %     t_ini=t+1;
    %     time = time
    
    if ~(exist('time','var')==1)
        time =t*model.advection.dt_adv;
    end
    
    day_last_plot = floor(time/24/3600);
else
    %     t_ini=1;
    
    
    if model.advection.cov_and_abs_diff
        cov_w = nan(1,N_t);
        day_ref_cov = 100;
        %day_ref_cov = 1;
        %warning(['This reference time for the covariance should ' ...
        %    'correspond to the stationnary regime']);
        % w_save = w;
        load( [model.folder.folder_simu '/files/' num2str(day_ref_cov)...
            '.mat'],'w','t');
        w_ref_cov = w(:)'; clear w;
        t_ref_cov = t; clear t;
        % w = w_save; clear w_save;
    end
    
    time = 0;
    w_fv = w;
    
    day_last_plot = -inf;
end

% if model.sigma.sto
%     for sampl=1:N_ech
%         model_sampl(sampl)=model;
%     end
% end

%%

while time < model.advection.advection_duration
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_b);
    w = real(ifft2(fft_w));
    % clear fft_w
    
    %% Comp ADSD Self Sim and EOF
    day_num = (floor(time/24/3600));
%     if (first_time == time) && ...
    if (day_num == 105) && plt_first_time && ...
            strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        plt_first_time = false;
        color_sexy(1,:) = [0.8500, 0.3250, 0.0980];
        color_sexy(2,:) = [0.9290, 0.6940, 0.1250];
        color_sexy(3,:) = [0.4940, 0.1840, 0.5560];
        color_sexy(4,:) = [0.4660, 0.6740, 0.1880] 	;
        color_sexy(5,:) = [0.3010, 0.7450, 0.9330];
        color_sexy(6,:) = [0.6350, 0.0780, 0.1840];
        LineStyle_ =[ '-' ];
        Marker_ = {'none','o','+','*','x','s'};
        
        %     plot_abs_diff_from_sigma_postprocess(model,fft2(sigma_dBt_on_sq_dt))
        fct_sigma_spectrum_abs_diff_postprocess(model,fft2(w),true,'0',...
            color_sexy(6,:));
        
        % Load precomputed EOFs and correspond variance tensor
        %     current_folder = pwd;
        %     cd(model.folder.folder_simu);
        %     cd ..
        %     model.folder.folder_EOF = [ pwd '/folder_EOF'];
        model.folder.folder_EOF = [ pwd '/images/SQG_MU_HV_4/' ...
            'type_spectrum_sigma_EOF/' ...
            'disym_Vortices_forced_turb_Spring/64x64/folder_EOF'];
        %     cd(current_folder);
        load([ model.folder.folder_EOF '/EOF.mat'],'EOF');
        EOF = permute(EOF,[1 2 3 5 4]);
        sigma = EOF; clear EOF;
        
        nb_EOF_v = [8000 2000 200 20 2];
        for k=1:length(nb_EOF_v)
            nb_EOF = nb_EOF_v(k);
            sigma = sigma(:,:,:,:,1:nb_EOF);
            %         if model.sigma.nb_EOF > 400
            %             sigma_dBt_on_sq_dt = 0;
            %             for j=1: nb_EOF
            %                 sigma_dBt_on_sq_dt = sigma_dBt_on_sq_dt + ...
            %                     sigma(:,:,:,:,j) .* ...
            %                     randn( [ 1 1 1 N_ech ]);
            %             end
            %         else
            %             sigma_dBt_on_sq_dt = sum( sigma .* ...
            %                 randn( [ 1 1 1 N_ech nb_EOF ]) , 5);
            %         end
            sigma_dBt_on_sq_dt = sum( sigma .* ...
                randn( [ 1 1 1 10 nb_EOF ]) , 5);
            plot_abs_diff_from_sigma_postprocess_add(model,...
                fft2(sigma_dBt_on_sq_dt), color_sexy(end-k,:), ...
                LineStyle_, Marker_{k+1} );
            %             fft2(sigma_dBt_on_sq_dt), [0.8 0.1 (k-1)/5]);
            %         %             fft2(sigma_dBt_on_sq_dt), [0.8 0.1 0.1 + (k-1)/5]);
        end
        pl = findobj(gca,'Type','line');
        %     subplot(1,2,1);
        figure10=figure(10);
        width = 9;
        height = 3;
        set(figure10,'Units','inches', ...
            'Position',[0 0  width height], ...
            'PaperPositionMode','auto');
        axP = get(gca,'Position');
        lgd=legend(pl([10 9 5 4 3 2 1 ]),...
            {'w',...
            '$\sigma \dot{B} \  $   Self.Sim.',...
            '$\sigma \dot{B} \ $   8000 EOFs',...
            '$\sigma \dot{B} \ $   2000 EOFs',...
            '$\sigma \dot{B} \ $   200 EOFs',...
            '$\sigma \dot{B} \ $   20 EOFs',...
            '$\sigma \dot{B} \ $   2 EOFs'},...
            'Interpreter','latex',...
            'Location','northwestoutside');
        %         'Location','northeastoutside');
        % %         'Location','northeastoutside');
        set(gca, 'Position', axP);
        set(gcf,'children',flipud(get(gcf,'children')));
        %     set(figure10,'Units','inches', ...
        %         'Position',[0 0  width*1.25 height], ...
        %         'PaperPositionMode','auto');
        v=sigma_dBt_on_sq_dt;
        2*mean(v(:).^2)
        drawnow
        
        folder_simu = model.folder.folder_simu;
        eval( ['print -depsc ' folder_simu ...
            '/Comp_ADSD_SelfSim_EOF.eps']);
        keyboard;
        %     v=sigma_dBt_on_sq_dt/sqrt(model.advection.dt_adv);
        %     2*mean(v(:).^2)*model.advection.dt_adv
        
        %    taille_police = 8;
        % %     legend('-7/3','fit','w','\sigma \dot{B}',[],[],[],...
        %     legend({'-7/3','fit','w','$\sigma \dot{B}$','','','',...
        %         '8000 EOF','2000 EOF','200 EOF','20 EOF','2 EOF'},...
        %         'FontUnits','points',...
        %         'FontWeight','normal',...
        %         'FontSize',taille_police,...
        %         'interpreter','latex',...
        %         'FontName','Times');
    end
    %%
    
    
    %% Time-uncorrelated velocity (isotropic and homogeneous in space)
    if ~ model.sigma.sto % Deterministic case
        sigma_dBt_on_sq_dt = 0;
    else % Stochastic case
        
        if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            %             sigma = nan([model.grid.MX 2 N_ech]);
            %             tr_a = nan(1,N_ech);
            %             slope_sigma = nan(1,N_ech);
            %             offset_spectrum_a_sigma = nan(1,N_ech);
            %             km_LS = nan(1,N_ech);
            %%
            % for sampl=1:N_ech
            % parfor sampl=1:N_ech
            %%
            
%             if model.sigma.time_smooth.bool
%                 sigma_s = sigma;
%             end
            
            % [sigma(:,:,:,sampl), ~, tr_a(sampl) ] ...
            [ sigma, ~, model.sigma.tr_a ,....
                model.sigma.slope_sigma,...
                model.sigma.offset_spectrum_a_sigma, ...
                model.sigma.km_LS ]...
                = fct_sigma_spectrum_abs_diff_mat(...
                model,fft_w,false);
            
            if model.sigma.time_smooth.bool
                sigma_n_s = sigma;
                d_sigma_s = 1/model.sigma.time_smooth.tau * ...
                                 ( - sigma_s + sigma_n_s ) ;
                sigma_s = sigma_s + d_sigma_s * model.advection.dt_adv;
                sigma = sigma_s;
            end
            
            %             %                 % [sigma(:,:,:,sampl), ~, tr_a(sampl) ] ...
            %             %                 [ sigma(:,:,:,sampl), ~, tr_a(sampl) ,....
            %             %                     slope_sigma(sampl),...
            %             %                     offset_spectrum_a_sigma(sampl), ...
            %             %                     km_LS(sampl) ]...
            %             %                     = fct_sigma_spectrum_abs_diff(...
            %             %                     model,fft_w(:,:,:,sampl),false);
            model.sigma.a0 = model.sigma.tr_a/2;
            model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
            % Diffusion coefficient
            model.advection.coef_diff = permute( 1/2 * model.sigma.a0 , ...
                [1 3 4 2]);
            if model.sigma.assoc_diff | model.sigma.Smag.bool
                % warning('deal with slope when there are several realizations')
                model.sigma.a0_SS = ...
                    1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
                    model, ...
                    model.sigma.kappaMinUnresolved_on_kappaShanon ...
                    *(pi/sqrt(prod(model.grid.dX))), ...
                    model.sigma.kappaMaxUnresolved_on_kappaShanon ...
                    *(pi/sqrt(prod(model.grid.dX))));
                model.sigma.a0_LS = ...
                    model.sigma.a0 ;
                model.sigma.a0 = ...
                    model.sigma.a0 ...
                    + model.sigma.a0_SS;
                model.sigma.a0_on_dt = model.sigma.a0 ...
                    / model.advection.dt_adv;
                % Diffusion coefficient
                model.advection.coef_diff = permute(...
                    1/2 * model.sigma.a0 ,[1 3 4 2]);
                %end
                
                
                if model.sigma.Smag.bool
                    iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
                    %if model.sigma.a0_SS > eps
                    if any(iii_non_degenerate_a0_SS)
                        if model.sigma.Smag.epsi_without_noise
                            sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                                bsxfun(@times,...
                                permute( ...
                                sqrt(2./(...
                                model.sigma.a0_LS(iii_non_degenerate_a0_SS) ...
                                + model.sigma.a0_SS(iii_non_degenerate_a0_SS))) ...
                                , [ 1 3 4 2]) , ...
                                sigma(:,:,:,(iii_non_degenerate_a0_SS)) );
                        else
                            sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                                bsxfun( @times, ...
                                permute( ...
                                sqrt(2./...
                                model.sigma.a0_SS(iii_non_degenerate_a0_SS)) ...
                                , [ 1 3 4 2]) , ...
                                sigma(:,:,:,iii_non_degenerate_a0_SS) );
                        end
                        %                     else
                        %                         sigma = [];
                    end
                    if any(~iii_non_degenerate_a0_SS)
                        if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                            % The absolute diffusivity diagnosed from the large-scale
                            % kinematic spectrum is too weak. It suggests that there
                            % are few small scales and no subgrid terms is needed.
                            % Moreover, setting subgris terms to zero prevent numerical
                            % errors.
                            sigma(:,:,:,~iii_non_degenerate_a0_SS) = ...
                                zeros(size(sigma(:,:,:,~iii_non_degenerate_a0_SS)));
                            model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
                        else
                            error('Unknow case');
                        end
                    end
                    %                     if model.sigma.a0_SS > eps
                    %                         if model_sampl(sampl).sigma.Smag.epsi_without_noise
                    %                             sigma(:,:,:,sampl) = ...
                    %                                 sqrt(2/(model_sampl(sampl).sigma.a0_LS ...
                    %                                 + model_sampl(sampl).sigma.a0_SS)) ...
                    %                                 * sigma(:,:,:,sampl);
                    %                         else
                    %                             sigma(:,:,:,sampl) = ...
                    %                                 sqrt(2/model_sampl(sampl).sigma.a0_SS) ...
                    %                                 * sigma(:,:,:,sampl);
                    %                         end
                    %                     elseif strcmp(model_sampl(sampl).sigma.type_spectrum,'SelfSim_from_LS')
                    %                         % The absolute diffusivity diagnosed from the large-scale
                    %                         % kinematic spectrum is too weak. It suggests that there
                    %                         % are few small scales and no subgrid terms is needed.
                    %                         % Moreover, setting subgris terms to zero prevent numerical
                    %                         % errors.
                    %                         sigma(:,:,:,sampl) = ...
                    %                             zeros(size(sigma(:,:,:,sampl)));
                    %                         model_sampl(sampl).advection.coef_diff = 0;
                    %                     else
                    %                         error('Unknow case');
                    %                     end
                end
                
            end
        elseif model.sigma.Smag.bool | model.sigma.assoc_diff
            model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
        else
            % Variance tensor
            %parfor sampl=1:N_ech
            model.sigma.a0 = 2 * model.physical_constant.f0 ...
                / model.sigma.k_c^2;
            %end
        end
        
        
        if model.sigma.hetero_energy_flux
            
            %parfor sampl=1:N_ech
            % for sampl=1:N_ech
            
            coef_modulation = ...
                fct_epsilon_k_onLine(model,fft_b,fft_w);
            %end
        elseif model.sigma.hetero_modulation | ...
                model.sigma.hetero_modulation_V2
            if ( isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                    &&    model.sigma.hetero_energy_flux_prefilter ) ...
                    || (isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                    &&    model.sigma.hetero_energy_flux_postfilter)
                error('not coded yet')
            end
            coef_modulation = ...
                fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
        elseif model.sigma.hetero_modulation_Smag
            % Heterogeneous dissipation coefficient
            if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                    &&    model.sigma.hetero_energy_flux_prefilter
                % Pre-filtering
                fft_b_for_modulation = bsxfun(@times, ...
                    model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                    fft_b);
            else
                fft_b_for_modulation = fft_b;                
            end
            coef_modulation = fct_coef_diff(model,fft_b_for_modulation);
            clear fft_b_for_modulation
            
            if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                    &&    model.sigma.hetero_energy_flux_postfilter
                % Post-filtering
                coef_modulation = real(ifft2(bsxfun(@times, ...
                    model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                    fft2(coef_modulation))));             
            end
            m_coef_modulation = mean(mean(coef_modulation,1),2);
            coef_modulation = bsxfun( @times, ...
                1./m_coef_modulation, coef_modulation);
            clear m_coef_modulation
        elseif model.sigma.Smag.bool
            % Heterogeneous dissipation coefficient
            coef_modulation = fct_coef_diff(model,fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            coef_modulation = model.sigma.Smag.coef_Smag * coef_modulation ;
            
            %             %     figure(12);fct_spectrum( model,fft2(sigma_dBt_dt));
            %             %   figure(13);fct_spectrum( model,fft2(coef_diff_aa));
            %             %     figure(14);fct_spectrum( model,fft2(sqrt(coef_diff_aa)));
            %             %    figure(15);imagesc(sqrt(coef_diff_aa)');axis xy;axis equal
            %             %%
            %             % for sampl=1:N_ech
            %parfor sampl=1:N_ech
            
            % Coefficient coef_Smag to target a specific diffusive scale
            iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
            %if model_sampl(sampl).sigma.a0_SS > eps
            model.advection.coef_diff = zeros([1 1 1 model.advection.N_ech]);
            if model.sigma.Smag.epsi_without_noise
                model.advection.coef_diff(iii_non_degenerate_a0_SS) = 1;
            else
                % Taking into account the noise in the energy budget
                %                 model.advection.coef_diff(:,:,:,iii_non_degenerate_a0_SS) = ...
                model.advection.coef_diff(iii_non_degenerate_a0_SS) = ...
                    permute( ...
                    (1 + model.sigma.a0_LS(iii_non_degenerate_a0_SS) ./ ...
                    model.sigma.a0_SS(iii_non_degenerate_a0_SS)) , ...
                    [1 4 3 2] );
            end
            if any(~iii_non_degenerate_a0_SS)
                if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                    % The absolute diffusivity diagnosed from the large-scale
                    % kinematic spectrum is too weak. It suggests that there
                    % are few small scales and no subgrid terms is needed.
                    % Moreover, setting subgris terms to zero prevent numerical
                    % errors.
                    model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
                else
                    error('Unknow case');
                end
            end
%             %                 if model_sampl(sampl).sigma.a0_SS > eps
%             %                     if model_sampl(sampl).sigma.Smag.epsi_without_noise
%             %                         model_sampl(sampl).advection.coef_diff = 1;
%             %                     else
%             %                         % Taking into account the noise in the energy budget
%             %                         model_sampl(sampl).advection.coef_diff = ...
%             %                             (1 + model_sampl(sampl).sigma.a0_LS / ...
%             %                             model_sampl(sampl).sigma.a0_SS) ;
%             %                     end
%             %                 elseif strcmp(model_sampl(sampl).sigma.type_spectrum,'SelfSim_from_LS')
%             %                     % The absolute diffusivity diagnosed from the large-scale
%             %                     % kinematic spectrum is too weak. It suggests that there
%             %                     % are few small scales and no subgrid terms is needed.
%             %                     % Moreover, setting subgris terms to zero prevent numerical
%             %                     % errors.
%             %                     model_sampl(sampl).advection.coef_diff = 0;
%             %                 else
%             %                     error('Unknow case');
%             %                 end
%             model.advection.coef_diff = ...
%                 bsxfun( @times, ...
%                 model.advection.coef_diff,...
%                 coef_modulation) ;
%             %end
            
            if model.sigma.Smag.SS_vel_homo
                coef_modulation = mean(mean(coef_modulation,2),1);
                %coef_modulation = mean(coef_modulation(:));
            end
        else
            coef_modulation = 1;
        end
        %% Variance tensor
        model.advection.coef_diff = ...
            bsxfun( @times, ...
            model.advection.coef_diff,...
            coef_modulation) ;
        %%
        if model.sigma.assoc_diff
            
            % for sampl=1:N_ech
            %parfor sampl=1:N_ech
            model.sigma.a0 = ...
                model.sigma.a0_LS ...
                + model.sigma.a0_SS;
            model.sigma.a0_on_dt = ...
                model.sigma.a0 / ...
                model.advection.dt_adv;
            % Diffusion coefficient
            %                 model.advection.coef_diff = ...
            %                     coef_modulation * ...
            %                     1/2 * model.sigma.a0;
            model.advection.coef_diff = bsxfun(@times, ...
                permute( 1/2 * model.sigma.a0 , [ 1 3 4 2]) , ...
                coef_modulation ) ;
            %end
        end
        
        
        
        
        % Maximum of the variance tensor
        %         % coef_modulation_a0 = max(max(coef_modulation,[],2),[],1);
        %         a0_temp = nan([1 1 1 N_ech]);
        a0_temp = nan([N_ech 1]);
        if size(coef_modulation,4)==1
            coef_modulation = repmat(coef_modulation,[1 1 1 N_ech]);
        end
        
        % for sampl=1:N_ech
        %parfor sampl=1:N_ech
        
        model.sigma.a0 = bsxfun(@times, ...
            permute(model.sigma.a0 , [ 1 3 4 2]) , ...
            coef_modulation ) ;
        
        
        
        if model.sigma.Smag.bool
            iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
            %if model_sampl(sampl).sigma.a0_SS > eps
            if model.sigma.Smag.epsi_without_noise
                model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                    model.sigma.a0(iii_non_degenerate_a0_SS) ...
                    ./ (model.sigma.a0_SS(iii_non_degenerate_a0_SS) ...
                    + model.sigma.a0_LS(iii_non_degenerate_a0_SS));
            else
                % Taking into account the noise in the energy budget
                model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                    model.sigma.a0(iii_non_degenerate_a0_SS) ...
                    ./ model.sigma.a0_SS(iii_non_degenerate_a0_SS);
            end
            if any(~iii_non_degenerate_a0_SS)
                if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                    % The absolute diffusivity diagnosed from the large-scale
                    % kinematic spectrum is too weak. It suggests that there
                    % are few small scales and no subgrid terms is needed.
                    % Moreover, setting subgris terms to zero prevent numerical
                    % errors.
                    model.sigma.a0(~iii_non_degenerate_a0_SS) = 0;
                else
                    error('Unknow case');
                end
            end
        end
        
        a0_temp = max(max( model.sigma.a0 ,[],2), [],1);
        model.sigma.a0 = max(a0_temp(:)) ; clear a0_temp;
        
        %% Simulation of sigma dBt
        
        if ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...             
                strcmp(model.sigma.type_spectrum,'Euler_EOF'))
            % Fourier transform of white noise
            dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
            % Multiplication by the Fourier transform of the kernel \tilde \sigma
            fft_sigma_dBt_on_sq_dt = bsxfun(@times,sigma,dBt_C_on_sq_dt);
            clear dBt_C_on_sq_dt
            % Homogeneous velocity field
            sigma_dBt_on_sq_dt = real(ifft2(fft_sigma_dBt_on_sq_dt));
            clear fft_sigma_dBt_on_sq_dt
        else
            sigma_dBt_on_sq_dt = sum( sigma .* ...
                randn( [ 1 1 1 N_ech model.sigma.nb_EOF ]) , 5);
        end
        
        %         % Fourier transform of white noise
        %         dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
        %         % Multiplication by the Fourier transform of the kernel \tilde \sigma
        %         fft_sigma_dBt_on_dt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
        %         clear dBt_C_on_sq_dt
        %         % Homogeneous velocity field
        %         sigma_dBt_dt = real(ifft2(fft_sigma_dBt_on_dt));
        %         clear fft_sigma_dBt_on_dt
        
        % Heterogeneous small-scale velocity
        sigma_dBt_on_sq_dt = bsxfun(@times, sqrt(coef_modulation) , ...
            sigma_dBt_on_sq_dt);
        
        if model.sigma.proj_free_div & any(coef_modulation(:) ~= 1)
            % nrj_before_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
            sigma_dBt_on_sq_dt = fct_proj_free_div(model,sigma_dBt_on_sq_dt);
            % nrj_after_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
        end
        
    end
    
    %% Specify determinstic heterogeneous subgrid model
    if model.advection.Smag.bool
        if model.advection.Lap_visco.bool
            % Heterogeneous dissipation coefficient
            model.advection.coef_diff = fct_coef_diff(model,fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            model.advection.coef_diff = ...
                model.advection.Smag.coef_Smag * ...
                model.advection.coef_diff ;
        elseif model.advection.HV.bool
            % Heterogeneous HV coefficient
            model.advection.coef_diff = ...
                fct_coef_HV(model,fft_b);
        end
        % Maximum dissipation coefficient
        model.advection.HV.maxVal = max(model.advection.coef_diff(:));
    end
    
    %% Dynamical time step
    model.advection.dt_adv = fct_CFL(model,w);
    %     if model.sigma.sto
    %         %         model_sampl.advection.dt_adv = repamt( model.advection.dt_adv,...
    %         %             [1 1 1 N_ech]);
    %         % model_sampl.advection.dt_adv = model.advection.dt_adv;
    %         %%
    %         parfor sampl=1:N_ech
    %             % for sampl=1:N_ech
    %             %%
    %             model_sampl(sampl).advection.dt_adv = model.advection.dt_adv;
    %         end
    %     end
    time = time + model.advection.dt_adv;
    % for t=t_ini:N_t
    % for t=1:N_t
    
    %% Adding time-correlated and time decorrelated velocity
    % w_fv = w;
    
    
    if  model.sigma.sto & ~model.sigma.no_noise
        w = w + sigma_dBt_on_sq_dt/sqrt(model.advection.dt_adv);
    end
    % if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %         w(:,:,1) = w(:,:,1) + Vy;
    %     end
    
    %% Transport of tracer
    if ~ model.sigma.sto
        %if ~ model.sigma.sto
        % Runge-Kutta 4 scheme
        fft_b = RK4_fft_advection(model,fft_b, w);
    else
        % Euler scheme
        fft_b = fft_b ...
            + deriv_fft_advection( ...
            model, fft_b, w) ...
            * model.advection.dt_adv;
        clear model_temp
    end
    
    %% Discard particles which have blown up
    iii = isnan(fft_b) | isinf(abs(fft_b));
    if any(iii(:))
        iii=any(any(any(iii,3),2),1);
        if all(iii(:))
            if model.plots
                error('The simulation has blown up');
            else
                warning('One simulation has blown up');
                fprintf('Folder of the simulation : \n');
                fprintf([ model.folder.folder_simu ' \n']);
                return;
            end
        end
        nb_dead_pcl = sum(iii);
        warning([ num2str(nb_dead_pcl) ' particle(s) on ' num2str(N_ech) ...
            ' have(s) blown up and' ...
            ' are(is) resampled uniformly on' ...
            ' the set of the others particles']);
        N_ech_temp = N_ech - nb_dead_pcl;
        fft_b(:,:,:,iii)=[];
        iii_sample = randi(N_ech_temp,nb_dead_pcl,1);
        for k=1:nb_dead_pcl
            fft_b(:,:,:,end+1) = fft_b(:,:,:,iii_sample(k));
        end
    end
    clear iii
    
    %% Veloicity covariance and Eulerian absolute diffusivity
    if model.advection.cov_and_abs_diff
        cov_w(t) = 1/prod(model.grid.MX) * w_ref_cov * w(:) ;
    else
        cov_w = nan ;
    end
    
    %% Plots and save
    %     % tt = floor(t ); % Number of days
    %     tt = floor(t *model.advection.dt_adv/ (3600*24)); % Number of days
    %     tt = floor(time/ (3600*24)); % Number of days
    %     if tt > tt_last
    %         tt_last = tt;
    %         day = num2str(floor(time/24/3600));
    %         % day = num2str(floor(t*model.advection.dt_adv/24/3600));
    %         % t_last_plot = t;
    
    day_num = (floor(time/24/3600));
    
    %%
    %     warning('DEBUG')
    %     % if (t_loop - t_last_plot)*dt >= 3600*24*1
    %     day_num = (floor(time/24/3600));
    %     day = num2str(day_num);
    %     day
    %     day_last_plot = day_num;
    %
    %     nb=2;
    %     fft_w=fft_w(:,:,:,1:nb);
    %     model.advection.N_ech =nb;
    %     fct_sigma_spectrum_abs_diff_mat(model,fft_w,true,day);
    %%
    
    
    if day_num > day_last_plot
        % if (t_loop - t_last_plot)*dt >= 3600*24*1
        day_num = (floor(time/24/3600));
        day = num2str(day_num);
        day
        day_last_plot = day_num;
        
        if model.plots
            toc
            tic
            %             if ~model.sigma.sto
            %                 model_sampl = model;
            %             end
            
            fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
            a_0_LS = mean(sigma_dBt_on_sq_dt(:).^2);
            if model.sigma.sto
                a_0_LS
            end
            
            id_part=1;
            %%
            if model.advection.Smag.bool || ...
                    (model.sigma.sto && ( model.sigma.Smag.bool ...
                    || model.sigma.hetero_modulation ...
                    || model.sigma.hetero_modulation_V2 ...
                    || model.sigma.hetero_energy_flux ...
                    || model.sigma.hetero_modulation_Smag ) )
                % Coefficient coef_Smag to target a specific diffusive scale
                if model.advection.Smag.bool
                    [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                        fft_b(:,:,:,id_part));
                    coef_diff_aa = ...
                        model.advection.Smag.coef_Smag * coef_diff_aa ;
                    coef_diff = ...
                        model.advection.Smag.coef_Smag * coef_diff ;
                elseif model.sigma.hetero_modulation_Smag
                    if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                            &&    model.sigma.hetero_energy_flux_prefilter
                        % Pre-filtering
                        fft_b_for_modulation = bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                            fft_b);
                    else
                        fft_b_for_modulation = fft_b;
                    end
                    [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                        fft_b_for_modulation(:,:,:,id_part));
                    
                    
                    if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                            &&    model.sigma.hetero_energy_flux_postfilter
                        % Post-filtering
                        coef_modulation = real(ifft2(bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                            fft2(coef_modulation))));
                        coef_diff_aa = real(ifft2(bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                            fft2(coef_diff_aa))));
                    end
            
                    coef_diff_aa = coef_diff_aa / mean(coef_diff_aa(:));
                    coef_diff = coef_diff / mean(coef_diff(:));
                elseif model.sigma.Smag.bool
                    [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                        fft_b(:,:,:,id_part));
                    if model.sigma.a0_SS(id_part) > eps
                        coef_diff = ...
                            (1 + model.sigma.a0_LS(id_part) ...
                            / model.sigma.a0_SS(id_part)) * ...
                            model.sigma.Smag.coef_Smag * coef_diff ;
                        coef_diff_aa = ...
                            (1 + model.sigma.a0_LS(id_part) ...
                            / model.sigma.a0_SS(id_part)) * ...
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
                    if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                            &&    model.sigma.hetero_energy_flux_prefilter
                        % Pre-filtering
                        fft_b_for_modulation = bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                            fft_b);
                    else
                        fft_b_for_modulation = fft_b;
                    end
                    
                    fft_w_for_modulation = SQG_large_UQ(model, fft_b_for_modulation);
                    %                     coef_diff_aa = ...
                    %                         model.sigma.a0(id_part)/2  ;
                    [coef_diff_aa,coef_diff] = ...
                        fct_coef_estim_AbsDiff_heterogeneous(...
                        model,fft_w_for_modulation(:,:,:,id_part));
                    
                    if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                            &&    model.sigma.hetero_energy_flux_postfilter
                        % Post-filtering
                        coef_modulation = real(ifft2(bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                            fft2(coef_modulation))));
                        coef_diff_aa = real(ifft2(bsxfun(@times, ...
                            model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                            fft2(coef_diff_aa))));
                    end
            
                    coef_diff_aa = ...
                        model.sigma.a0(id_part)/2 * coef_diff_aa ;
                    coef_diff = ...
                        model.sigma.a0(id_part)/2 * coef_diff ;
                elseif model.sigma.hetero_energy_flux
                    %                     coef_diff_aa = ...
                    %                         model.sigma.a0(id_part)/2  ;
                    coef_diff_aa = ...
                        fct_epsilon_k_onLine(model,fft_b,fft_w);
                    coef_diff_aa = coef_diff_aa(:,:,:,id_part);
                    coef_diff_aa = ...
                        model.sigma.a0(id_part)/2 * coef_diff_aa ;
                    coef_diff=nan(size(coef_diff_aa));
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
                %%
                figure(9);
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
                fct_plot(model,fft_b,day);
            
            if model.advection.plot_dissip
                fct_plot_dissipation(model,fft_b,sigma,day);
            end
            
            if model.sigma.sto & ...
                    strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                fct_sigma_spectrum_abs_diff_mat(model,fft_w,true,day);
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
                slope_sigma=model.sigma.slope_sigma(id_part)
                a0_LS=model.sigma.a0_LS(id_part)
                a0_SS=model.sigma.a0_SS(id_part)
            end
            
            if model.advection.cov_and_abs_diff
                abs_diff = sum(cov_w(t_ref_cov:end))*model.advection.dt_adv;
                figure(36)
                plot(model.advection.dt_adv*(0:(N_t-1))/3600/24,cov_w);
                hold on;
                plot(t_ref_cov*[1 1]*model.advection.dt_adv/3600/24,...
                    max(abs(cov_w(~isnan(cov_w))))*[-1 1],'r');
                hold off
            end
            
            % Dissipation by scale
            if model.advection.plot_epsilon_k
                fct_plot_epsilon_k(model,fft_b(:,:,:,id_part),day);
                % fct_plot_epsilon_k(model,fft_b,int_epsilon,day);
            end
            dt = model.advection.dt_adv
            
        end
        
        
        % Save files
        save( [model.folder.folder_simu '/files/' day '.mat'], ...
            'model','time','fft_b','w','sigma_dBt_on_sq_dt', ...
            'sigma');
        %             'sigma_on_sq_dt');
        %         %             'model','t','fft_b','w','sigma_dBt_on_sq_dt', ...
        %         %             'sigma_on_sq_dt');
        %         %         %             'sigma_on_sq_dt','cov_w','abs_diff');
    end
    
    %     % Dissipation by scale
    %     if (t == t_last_plot + 1) &  model.advection.plot_epsilon_k
    %         fct_plot_epsilon_k(model,fft_b,day);
    %         % fct_plot_epsilon_k(model,fft_b,int_epsilon,day);
    %     end
    
end
toc