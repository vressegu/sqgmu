function [turb_dissip, noise_intake, estim_noise_intake, ...
    estim_aa_noise_intake,...
    dissip_HV, intake_forcing,model] = ...
    fct_dissip(model,fft_b,sigma_on_sq_dt)
% Compute the deterministic and stochastic dissipation terms locally in
% space.
%

%warning('This function is only valid for homogeneous dissipation');



%% Grid of wave vectors
ikx = 1i*model.grid.k.kx; %"normal" grid
iky = 1i*model.grid.k.ky;
k2 = model.grid.k.k2;
ZM = model.grid.k.ZM; %index of modes ot be zero'ed out
ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
iky_aa = model.grid.k_aa.iky;
k2_aa = model.grid.k_aa.k2;
mask_aa = model.grid.k_aa.mask; %anti-aliasing mask

%% Buoyancy
b = real(ifft2(fft_b));

%% Gradient of b
gradb(:,:,1,:) = real(ifft2( ikx.*fft_b ));
gradb(:,:,2,:) = real(ifft2( iky.*fft_b ));

%% Gradient of b, anti-aliased
% in Fourier space, de-aliased, then in physical space.
gradb_aa(:,:,1,:) = real(ifft2( ikx_aa.*fft_b ));
gradb_aa(:,:,2,:) = real(ifft2( iky_aa.*fft_b ));

%% Stochastic terms
if ~ model.sigma.sto % Deterministic case
    noise_intake = zeros(size(b));
    estim_noise_intake = zeros(size(b));
    turb_dissip = zeros(size(b));
    estim_aa_noise_intake = zeros(size(b));
else % Stochastic case
    %% Advection term
    %if model.advection.N_ech <= 10
    N_ech_local = 100;
    % Fourier transform of white noise
    dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech_local]));
    % Multiplication by the Fourier transform of the kernel \tilde \sigma
    
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        dt = model.advection.dt_adv;
        % Time-correlated velocity
        fft_w = SQG_large_UQ(model, fft_b);
        [sigma, ~, tr_a ] ...
            = fct_sigma_spectrum_abs_diff( model,fft_w,false);
        a0 = tr_a/2;
        sigma_on_sq_dt = (1/sqrt(dt)) * sigma; clear sigma
        model.sigma.a0 = a0;
        model.sigma.a0_on_dt = model.sigma.a0 / dt;
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
        
        %%
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
                if model.sigma.a0_SS > eps
                    if model.sigma.Smag.epsi_without_noise
                        sigma_on_sq_dt = ...
                            sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                            * sigma_on_sq_dt;
                        %                             model.advection.coef_diff = 1;
                    else
                        sigma_on_sq_dt = sqrt(2/model.sigma.a0_SS) ...
                            * sigma_on_sq_dt;
                        %                             model.advection.coef_diff = 1 + ...
                        %                                 model.sigma.a0_LS / model.sigma.a0_SS ;
                    end
                elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                    % The absolute diffusivity diagnosed from the large-scale
                    % kinematic spectrum is too weak. It suggests that there
                    % are few small scales and no subgrid terms is needed.
                    % Moreover, setting subgris terms to zero prevent numerical
                    % errors.
                    sigma_on_sq_dt = zeros(size(sigma_on_sq_dt));
                    model.advection.coef_diff = 0;
                else
                    error('Unknow case');
                end
            end
            %%
        end
        
        %%
        
    elseif model.sigma.Smag.bool | model.sigma.assoc_diff
        model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
    else
        % Variance tensor
        model.sigma.a0 = 2 * model.physical_constant.f0 ...
            / model.sigma.k_c^2;
        
    end
    
    fft_sigma_dBt_dt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
    clear dBt_C_on_sq_dt
    
    % warning('need case self similar');
    
    sigma_dBt_dt = real(ifft2(fft_sigma_dBt_dt));
    if model.sigma.Smag.bool
        % Heterogeneous dissipation coefficient
        %coef_diff_aa = fct_coef_diff(model,fft_buoy_part);
        coef_diff_aa = fct_coef_diff(model,nan,gradb_aa);
        % Coefficient coef_Smag to target a specific diffusive scale
        coef_diff_aa = model.sigma.Smag.coef_Smag * coef_diff_aa ;
        
        if model.sigma.Smag.SS_vel_homo
            coef_modulation = mean(coef_diff_aa(:));
        else
            coef_modulation = coef_diff_aa;
        end
        
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
        
        
    elseif model.sigma.hetero_modulation | ...
            model.sigma.hetero_modulation_V2
        coef_modulation = ...
            fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
    elseif model.sigma.hetero_energy_flux
        coef_modulation = ...
            fct_epsilon_k_onLine(model,fft_b,fft_w);
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
    sigma_dBt_dt = bsxfun(@times, sqrt(coef_modulation) , ...
        sigma_dBt_dt);
    if model.sigma.proj_free_div
        sigma_dBt_dt = fct_proj_free_div(model,sigma_dBt_dt);
        %             for sampl=1:N_ech_local
        %                 sigma_dBt_dt(:,:,:,sampl) = fct_proj_free_div(...
        %                     model,sigma_dBt_dt(:,:,:,sampl));
        %             end
    end
    fft_sigma_dBt_dt = fft2(sigma_dBt_dt); clear sigma_dBt_dt
    
    % dealisasing of the velocity: FT, apply mask, iFT
    sigma_dBt_aa = ...
        real(ifft2( bsxfun(@times, fft_sigma_dBt_dt, mask_aa) ));
    clear fft_sigma_dBt_dt
    
    sigma_dBt_dtgradT2 = nan([model.grid.MX 1 N_ech_local]);
    for sampl=1:N_ech_local
        %parfor sampl=1:N_ech_local
        % Advective term in physical space
        sigma_dBt_dtgradT2(:,:,:,sampl) = ...
            sum(bsxfun(@times,sigma_dBt_aa(:,:,:,sampl),gradb_aa),3);
        % NB : in the stochastic case, w included both continuous and
        % non-continuous components of the velocity
        
    end
    % Remove aliasing
    %         sigma_dBt_dtgradT2(PX(1)+1,:,:,sampl)=0;
    %         sigma_dBt_dtgradT21(:,PX(2)+1,:,sampl)=0;
    sigma_dBt_dtgradT2=fft2(sigma_dBt_dtgradT2);
    sigma_dBt_dtgradT2(ZM(1),:,:,:)=0;
    sigma_dBt_dtgradT2(:,ZM(2),:,:)=0;
    sigma_dBt_dtgradT2=real(ifft2(sigma_dBt_dtgradT2));
    
    noise_intake = model.advection.dt_adv/2 ...
        * mean(sigma_dBt_dtgradT2.^2,4);
    
    %end
    
    %% Laplacian diffusion term (from the stochastic material derivative)
    if model.sigma.Smag.bool
        %         % Heterogeneous dissipation coefficient
        %         coef_diff_aa = fct_coef_diff(model,nan,gradb_aa);
        %         % Coefficient coef_Smag to target a specific diffusive scale
        %         coef_diff_aa = model.sigma.Smag.coef_Smag * coef_diff_aa ;
        %         % Taking into account the noise in the energy budget
        %         coef_diff_aa= model.advection.coef_diff * coef_diff_aa;
        
        if model.advection.Smag.spatial_scheme
            gradb_dissip = gradient_mat_2_per(b,model.grid.dX);
        else
            gradb_dissip = gradb;
        end
        
        % Dissipation
        turb_dissip_Smag = coef_diff_aa .* sum(gradb_dissip.^2,3);
        % turb_dissip = coef_diff_aa .* sum(gradb.^2,3);
        % Taking into account the noise in the energy budget
        
        if model.sigma.Smag.epsi_without_noise
            turb_dissip= turb_dissip_Smag;
        else
            turb_dissip= (1 + model.sigma.a0_LS / model.sigma.a0_SS) ...
                * turb_dissip_Smag;
        end
        
        % Estimation of the noise intake
        %         estim_noise_intake = ...
        %             1/(model.sigma.a0_SS/model.sigma.a0_LS + 1 ) ...
        %             * turb_dissip;
        if model.sigma.Smag.epsi_without_noise
            estim_noise_intake = model.sigma.a0_LS ...
                / (model.sigma.a0_LS + model.sigma.a0_SS) ...
                * turb_dissip_Smag;
        else
            estim_noise_intake = (model.sigma.a0_LS / model.sigma.a0_SS) ...
                * turb_dissip_Smag;
        end
        
        % Anti-Dissipation
        turb_dissip_estim_aa_Smag = coef_diff_aa .* sum(gradb_aa.^2,3);
        if model.sigma.Smag.epsi_without_noise
            turb_dissip_estim_aa = turb_dissip_estim_aa_Smag;
            % Estimation of the noise intake with antialiasing
            estim_aa_noise_intake = model.sigma.a0_LS ...
                / (model.sigma.a0_LS + model.sigma.a0_SS) ...
                * turb_dissip_estim_aa_Smag;
        else
            turb_dissip_estim_aa= (1 + model.sigma.a0_LS / model.sigma.a0_SS) ...
                * turb_dissip_estim_aa_Smag;
            % Estimation of the noise intake with antialiasing
            estim_aa_noise_intake = (model.sigma.a0_LS / model.sigma.a0_SS) ...
                * turb_dissip_estim_aa_Smag;
        end
        
    elseif model.sigma.hetero_modulation
        
        % Dissipation
        turb_dissip = coef_modulation .* sum(gradb.^2,3);
        % Taking into account the noise in the energy budget
        turb_dissip= model.sigma.a0/2 * turb_dissip;
        
        % Estimation of the noise intake
        estim_noise_intake = turb_dissip;
        
        % Anti-Dissipation
        turb_dissip_estim_aa = model.sigma.a0/2 * ...
            coef_modulation .* sum(gradb_aa.^2,3);
        % Estimation of the noise intake with antialiasing
        estim_aa_noise_intake = turb_dissip_estim_aa;
        
    else
        % Possibly aliased term !!!
        turb_dissip = model.advection.coef_diff * sum(gradb.^2,3);
        % turb_dissip = model.advection.coef_diff * sum(gradb_aa.^2,3);
        
        % % adv2 = - model.advection.coef_diff * k2 .* fft_b ;
        % % turb_dissip = - b .* real(ifft2(adv2));
        
        estim_noise_intake = turb_dissip;
        
        % Possibly aliased term !!!
        estim_aa_noise_intake = model.advection.coef_diff * ...
            sum(gradb_aa.^2,3);
    end
end

%% Deterministic subgrid model

% Specify determinstic heterogeneous subgrid model
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

if model.advection.Smag.bool
    %     % Stable k2
    %     k2=model.grid.k.k2;
    
    %% Heterogeneous HV or diffusivity/viscosity coefficient
    if model.advection.Lap_visco.bool
        %% Heterogeneous diffusivity/viscosity coefficient
        %         % square of the buoyancy gradient norm
        %         gradb2_aa = sum(gradb_aa.^2,3) ;
        %
        % %         % Root square of the buoyancy gradient norm
        % %         coef_diff = gradb2_aa .^(1/4) ;
        % %
        % %         % Coefficient coef_Smag to target a specific diffusive scale
        % %         coef_diff = model.advection.Smag.coef_Smag * coef_diff ;
        % %
        % %         % dealisasing of the diffusivity coefficient
        % %         coef_diff = fft2(coef_diff );
        % %         coef_diff_aa = real(ifft2( bsxfun(@times, coef_diff, mask_aa) ));
        % %         model.advection.HV.val = coef_diff_aa;
        %
        %
        %         % Heterogeneous dissipation coefficient
        %         coef_diff_aa = fct_coef_diff(model,fft_b,gradb_aa);
        %         % Coefficient coef_Smag to target a specific diffusive scale
        %         coef_diff_aa = model.advection.Smag.coef_Smag * coef_diff_aa ;
        %         model.advection.HV.val = coef_diff_aa;
        %
        %         % Dissipation
        %         dissip_HV = coef_diff_aa .* gradb2_aa;
        
        
        % Heterogeneous dissipation coefficient
        coef_diff_aa = fct_coef_diff(model,nan,gradb_aa);
        % Coefficient coef_Smag to target a specific diffusive scale
        coef_diff_aa  = model.advection.Smag.coef_Smag * coef_diff_aa  ;
        
        if model.advection.Smag.spatial_scheme
            gradb_dissip = gradient_mat_2_per(b,model.grid.dX);
        else
            gradb_dissip = gradb;
        end
        
        % Dissipation
        dissip_HV = coef_diff_aa .* sum(gradb_dissip.^2,3);
        
    elseif model.advection.HV.bool
        %% Heterogeneous HV coefficient
        % Laplacian at the power p of buoyancy
        Lap_p_b_aa = (-k2_aa) .^ (model.advection.HV.order/4) .* fft_b;
        Lap_p_b_aa = real(ifft2(Lap_p_b_aa));
        
        coef_HV = sqrt(abs( Lap_p_b_aa ));
        
        % Coefficient coef_Smag to target a specific diffusive scale
        coef_HV = model.advection.Smag.coef_Smag * coef_HV ;
        
        % dealisasing of the HV coefficient
        coef_HV = fft2(coef_HV );
        coef_HV_aa = real(ifft2( bsxfun(@times, coef_HV, mask_aa) ));
        model.advection.HV.val = coef_HV_aa;
        
        % Dissipation
        dissip_HV = coef_HV_aa .* Lap_p_b_aa .^2 ;
        
    else
        error('Unknown deterministic subgrid tensor');
    end
    
else
    %% Homogeneous hyper-viscosity
    % Possibly aliased term !!!
    
    %     % Unstable k2
    %     k2_HV=model.grid.k_HV.k2;
    
    % Laplacian at the power p of buoyancy
    adv4 = k2 .^ (model.advection.HV.order/4) .* fft_b;
    adv4 = real(ifft2(adv4));
    
    % Dissipation
    dissip_HV = model.advection.HV.val * adv4.^2;
    
end

%% Forcing
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            fft_forcing = model.advection.forcing.F;
        case 'Spring'
            fft_forcing =  ...
                - model.advection.forcing.on_T * ...
                ( fft_b - model.advection.forcing.F);
    end
    intake_forcing = b .* real(ifft2(fft_forcing));
else
    intake_forcing = 0;
end


