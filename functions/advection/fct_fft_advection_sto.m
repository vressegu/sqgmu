function [fft_b, model] = fct_fft_advection_sto(model,  fft_b)
% Advection of buoyancy using SQG or SQG_MU model
%

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
    elseif ( model.sigma.hetero_modulation |  ...
            model.sigma.hetero_modulation_V2 | ...
            model.sigma.hetero_energy_flux )
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
    if model.advection.N_ech >1
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
else
    % Variance tensor
    model.sigma.a0 = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
    % Diffusion coefficient
    model.advection.coef_diff = 1/2 * model.sigma.a0;
end


if model.sigma.sto & model.sigma.hetero_energy_flux
    coef_modulation = fct_epsilon_k_onLine(model,fft_b);
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
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %     Ly = model.grid.MX(2) * model.grid.dX(2);
    %     [~,Y]=ndgrid(model.grid.x,model.grid.y);
    %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         sin( 2 * pi * model.advection.forcing.freq_f * Y / Ly);
    %     clear Ly X Y
    
    
    model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
    model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
    [X,Y]=ndgrid(model.grid.x,model.grid.y);
    
    switch model.dynamics
        case 'SQG'
            U_caract = ...
                model.odg_b / model.physical_constant.buoyancy_freq_N;
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
                * model.advection.forcing.freq_f(1) * X ...
                + 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y );
            F = F / model.advection.forcing.on_T;
        case 'Spring'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing ...
                * sin( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X ) ...
                .* sin( 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y );
            %     F = model.advection.forcing.ampli_forcing * ...
            %         sin( 2 * pi * model.advection.forcing.freq_f * Y / ...
            %         model.advection.forcing.Ly);
    end
    
    
    
    model.advection.forcing.F = fft2(F);
    
    if strcmp(model.type_data, 'Zero')
        fft_w = SQG_large_UQ(model,  ...
            model.odg_b / model.advection.forcing.ampli_forcing ...
            * model.advection.forcing.F);
        %   w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
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
    clear Lx Ly X Y on_T U_caract ampli_scale
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
        a0 = tr_a/2;
        model.sigma.a0 = a0;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        missed_var_small_scale_spectrum = 2*a0;
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
    else
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
%                     sigma_on_sq_dt = ...
%                         sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
%                         * sigma_on_sq_dt;
                    sigma = ...
                        sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                        * sigma;
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
    else
        % Muliplicative constant of the kernel \tilde sigma
        sigma = ...
            sqrt(2*model.sigma.a0/missed_var_small_scale_spectrum) ...
            * sigma;
%         model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
%         sigma_on_sq_dt = ...
%             sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
%             * sigma; clear sigma
%         %             * sigma_on_sq_dt;
%         % the factor d=2 is for the dimension d of the space R^d
%         % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
        
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

%% Use a saved files of a former simulation ?
if model.advection.use_save
    warning(['The run begin from an older file instead of from the' ...
        'initial condition']);
    day = num2str(model.advection.day_save);
    name_file = [model.folder.folder_simu '/files/' day '.mat']; clear day
    load(name_file)
    
    % Version of matlab
    vers = version;
    year = str2double(vers(end-5:end-2));
    subvers = vers(end-1);
    model.folder.colormap_freeze = ...
        (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
    
    if ~isfield(model,'plots')
        model.plots = true;
    end
    
    t_ini=t+1;
    tt_last = -inf;
else
    t_ini=1;
    
    
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
    
    
end

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

time = 0;
w_fv = w;
while time < model.advection.advection_duration
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_b);
    w=real(ifft2(fft_w));
    % clear fft_w
    
    %% Time-uncorrelated velocity (isotropic and homogeneous in space)
    if ~ model.sigma.sto % Deterministic case
        sigma_dBt_dt = 0;
    else % Stochastic case
        
        if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            sigma = nan([model.grid.MX 2 N_ech]);
            tr_a = nan(1,N_ech);
            slope_sigma = nan(1,N_ech);
            offset_spectrum_a_sigma = nan(1,N_ech);
            km_LS = nan(1,N_ech);
            for sampl=1:N_ech
            % parfor sampl=1:N_ech
                % [sigma(:,:,:,sampl), ~, tr_a(sampl) ] ...
                [ sigma(:,:,:,sampl), ~, tr_a(sampl) ,....
                    slope_sigma(sampl),...
                    offset_spectrum_a_sigma(sampl), ...
                    km_LS(sampl) ]...
                    = fct_sigma_spectrum_abs_diff(...
                    model,fft_w(:,:,:,sampl),false);
            end
            a0 = tr_a/2;
            for sampl=1:N_ech
                model.sigma.slope_sigma(sampl) = slope_sigma(sampl);
                model.sigma.offset_spectrum_a_sigma(sampl) = ...
                    offset_spectrum_a_sigma(sampl);
                model.sigma.km_LS(sampl) = km_LS(sampl);
            end
            clear tr_a slope_sigma offset_spectrum_a_sigma km_LS
%             sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) ...
%                 * sigma; clear sigma
            model.sigma.a0 = a0;
            model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
            % Diffusion coefficient
            model.advection.coef_diff = 1/2 * model.sigma.a0;
            if model.sigma.assoc_diff | model.sigma.Smag.bool
                % warning('deal with slope when there are several realizations')
                model.sigma.a0_SS = ...
                    1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
                    model, ...
                    model.sigma.kappaMinUnresolved_on_kappaShanon ...
                    *(pi/sqrt(prod(model.grid.dX))), ...
                    model.sigma.kappaMaxUnresolved_on_kappaShanon ...
                    *(pi/sqrt(prod(model.grid.dX))));
                model.sigma.a0_LS = a0 ;
                model.sigma.a0 = a0 + model.sigma.a0_SS;
                model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
                % Diffusion coefficient
                model.advection.coef_diff = 1/2 * model.sigma.a0;
                
                %%
                if model.sigma.Smag.bool
%                     if model.sigma.a0_SS > eps
%                         if model.sigma.Smag.epsi_without_noise
%                             sigma_on_sq_dt = ...
%                                 sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
%                                 * sigma_on_sq_dt;
%                             %                             model.advection.coef_diff = 1;
%                         else
%                             sigma_on_sq_dt = sqrt(2/model.sigma.a0_SS) ...
%                                 * sigma_on_sq_dt;
%                             %                             model.advection.coef_diff = 1 + ...
%                             %                                 model.sigma.a0_LS / model.sigma.a0_SS ;
%                         end
%                     elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
%                         % The absolute diffusivity diagnosed from the large-scale
%                         % kinematic spectrum is too weak. It suggests that there
%                         % are few small scales and no subgrid terms is needed.
%                         % Moreover, setting subgris terms to zero prevent numerical
%                         % errors.
%                         sigma_on_sq_dt = zeros(size(sigma_on_sq_dt));
%                         model.advection.coef_diff = 0;
%                     else
%                         error('Unknow case');
%                     end
                    if model.sigma.a0_SS > eps
                        if model.sigma.Smag.epsi_without_noise
                            sigma = ...
                                sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                                * sigma;
                            %                             model.advection.coef_diff = 1;
                        else
                            sigma = sqrt(2/model.sigma.a0_SS) ...
                                * sigma;
                            %                             model.advection.coef_diff = 1 + ...
                            %                                 model.sigma.a0_LS / model.sigma.a0_SS ;
                        end
                    elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                        % The absolute diffusivity diagnosed from the large-scale
                        % kinematic spectrum is too weak. It suggests that there
                        % are few small scales and no subgrid terms is needed.
                        % Moreover, setting subgris terms to zero prevent numerical
                        % errors.
                        sigma = zeros(size(sigma));
                        model.advection.coef_diff = 0;
                    else
                        error('Unknow case');
                    end
                end
                %%
            end
        elseif model.sigma.Smag.bool | model.sigma.assoc_diff
            model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
        else
            % Variance tensor
            model.sigma.a0 = 2 * model.physical_constant.f0 ...
                / model.sigma.k_c^2;
        end
        
        
        % Fourier transform of white noise
        dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
        % Multiplication by the Fourier transform of the kernel \tilde \sigma
        fft_sigma_dBt_on_sq_dt = bsxfun(@times,sigma,dBt_C_on_sq_dt);
        clear dBt_C_on_sq_dt
        % Homogeneous velocity field
        sigma_dBt_on_sq_dt = real(ifft2(fft_sigma_dBt_on_sq_dt));
        clear fft_sigma_dBt
        
        if model.sigma.hetero_energy_flux
            coef_modulation = ...
                fct_epsilon_k_onLine(model,fft_b,fft_w);
        elseif model.sigma.hetero_modulation | ...
                model.sigma.hetero_modulation_V2
            coef_modulation = ...
                fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
        elseif model.sigma.Smag.bool
            % Heterogeneous dissipation coefficient
            coef_modulation = fct_coef_diff(model,fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            coef_modulation = model.sigma.Smag.coef_Smag * coef_modulation ;
            
            %     figure(12);fct_spectrum( model,fft2(sigma_dBt_on_sq_dt));
            %   figure(13);fct_spectrum( model,fft2(coef_diff_aa));
            %     figure(14);fct_spectrum( model,fft2(sqrt(coef_diff_aa)));
            %    figure(15);imagesc(sqrt(coef_diff_aa)');axis xy;axis equal
            
            % Coefficient coef_Smag to target a specific diffusive scale
            if model.sigma.a0_SS > eps
                if model.sigma.Smag.epsi_without_noise
                    model.advection.coef_diff = 1;
                else
                    % Taking into account the noise in the energy budget
                    model.advection.coef_diff = ...
                        (1 + model.sigma.a0_LS / model.sigma.a0_SS) ;
                end
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.advection.coef_diff = 0;
            else
                error('Unknow case');
            end
            model.advection.coef_diff = ...
                model.advection.coef_diff * coef_modulation;
            
            if model.sigma.Smag.SS_vel_homo
                coef_modulation = mean(mean(coef_modulation,2),1);
                %coef_modulation = mean(coef_modulation(:));
            end
        else
            coef_modulation = 1;
        end
        
        if model.sigma.assoc_diff
            %             model.sigma.a0_SS = coef_modulation * model.sigma.a0_SS;
            %             model.sigma.a0_LS = coef_modulation * model.sigma.a0_LS ;
            model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
            model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
            % Diffusion coefficient
            model.advection.coef_diff = coef_modulation * ...
                1/2 * model.sigma.a0;
            % model.advection.coef_diff = 1/2 * model.sigma.a0;
        end
        
        % Heterogeneous small-scale velocity
        sigma_dBt_on_sq_dt = bsxfun(@times, sqrt(coef_modulation) , ...
            sigma_dBt_on_sq_dt);
        
        if model.sigma.proj_free_div
            % nrj_before_proj_div = mean(sigma_dBt_dt(:).^2)
            sigma_dBt_on_sq_dt = fct_proj_free_div(model,sigma_dBt_on_sq_dt);
            % nrj_after_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
        end
        
        % Maximum of the variance tensor
        model.sigma.a0 = bsxfun( @times, model.sigma.a0 , coef_modulation );
        %%
        if model.sigma.Smag.bool
            if model.sigma.a0_SS > eps
                if model.sigma.Smag.epsi_without_noise
                    model.sigma.a0 = model.sigma.a0 ...
                        / (model.sigma.a0_SS + model.sigma.a0_LS);
                else
                    % Taking into account the noise in the energy budget
                    model.sigma.a0 = model.sigma.a0 / model.sigma.a0_SS;
                end
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.sigma.a0 = 0;
            else
                error('Unknow case');
            end
            
            %%
        end
        model.sigma.a0 = max(model.sigma.a0(:)) ;
        %model.sigma.a0 = model.sigma.a0 * max(coef_modulation(:))
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
        model_temp = cell(N_ech,1);
        for sampl=1:N_ech
        % parfor sampl=1:N_ech
            % Euler scheme
            model_temp{sampl} = model;
            if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                model_temp{sampl}.sigma.a0 = model.sigma.a0(sampl);
                % model_temp{sampl}.sigma.a0_on_dt =model.sigma.a0_on_dt(sampl);
                model_temp{sampl}.advection.coef_diff = ...
                   model.advection.coef_diff(sampl);
%                 model_temp{sampl}.advection.coef_diff = ...
%                     model.sigma.a0(sampl)/2 ...
%                     * coef_modulation(:,:,:,sampl) ;
            end
            fft_b(:,:,:,sampl) = fft_b(:,:,:,sampl) ...
                + deriv_fft_advection( ...
                model_temp{sampl}, fft_b(:,:,:,sampl), w(:,:,:,sampl)) ...
                * model.advection.dt_adv;
            %model, fft_b(:,:,:,sampl), w(:,:,:,sampl)) ...
            %   * model.advection.dt_adv;
            model_temp{sampl} = {};
        end
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
    tt = floor(time/ (3600*24)); % Number of days
    if tt > tt_last
        tt_last = tt;
        day = num2str(floor(time/24/3600));
        % day = num2str(floor(t*model.advection.dt_adv/24/3600));
        % t_last_plot = t;
        if model.plots
            fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
            if model.sigma.sto
                a_0_LS = mean(sigma_dBt_on_sq_dt(:).^2);
                a_0_LS
            end
            
            %%
            if model.advection.Smag.bool || ...
                    (model.sigma.sto && ( model.sigma.Smag.bool ...
                    || model.sigma.hetero_modulation ...
                    || model.sigma.hetero_modulation_V2 ) )
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
                fct_plot(model,fft_b,day);
            
            if ~model.sigma.sto
                sigma = 0;
                sigma_dBt_on_sq_dt=0;
            end
            if model.advection.plot_dissip
                fct_plot_dissipation(model,fft_b,sigma,day);
            end
            
            if model.sigma.sto & ...
                    strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                fct_sigma_spectrum_abs_diff(model,fft_w,true,day);
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
                fct_plot_epsilon_k(model,fft_b,day);
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