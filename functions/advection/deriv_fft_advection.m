function d_fft_b_adv = deriv_fft_advection(model, fft_b, w)
% Compute the Fourier transform of the partial derivative along time
% of the tracer according to the (stochastic or deterministic) transport
% equation with a hyperviscosity term
%

% Used of the branch of Pierre Derian (sqgmu2) for the de-aliasing

%% Grid of wave vectors
% kx = model.grid.k.kx; %"normal" grid
% ky = model.grid.k.ky;
k2 = model.grid.k.k2;
ZM = model.grid.k.ZM; %index of modes ot be zero'ed out
ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
iky_aa = model.grid.k_aa.iky;
k2_aa = model.grid.k_aa.k2;
mask_aa = model.grid.k_aa.mask; %anti-aliasing mask
% if model.advection.Smag.bool
% %     ikx_aa_LS = model.grid.k_aa_LS.ikx; %anti-aliased grid for gradient
% %     iky_aa_LS = model.grid.k_aa_LS.iky;
% %     k2_aa_LS = model.grid.k_aa_LS.k2;
%     mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask
% end
ikx = 1i*model.grid.k.kx;
iky = 1i*model.grid.k.ky;

%% Gradient of b
gradb(:,:,1,:) = real(ifft2( bsxfun(@times, ikx, fft_b )));
gradb(:,:,2,:) = real(ifft2( bsxfun(@times, iky, fft_b )));

%% Gradient of b, anti-aliased
% in Fourier space, de-aliased, then in physical space.
gradb_aa(:,:,1,:) = real(ifft2( bsxfun(@times, ikx_aa, fft_b )));
gradb_aa(:,:,2,:) = real(ifft2( bsxfun(@times, iky_aa, fft_b )));

%% Advection term

% % Gradient in Fourier space
% adv1x = 1i * kx .* fft_b;
% adv1y = 1i * ky .* fft_b;
%
% % Gradient in physical space
% gradb(:,:,1,:)=real(ifft2(adv1x));
% gradb(:,:,2,:)=real(ifft2(adv1y));

% dealisasing of the velocity: FT, apply mask, iFT
ft_w = fft2( w );
w_aa = real(ifft2( bsxfun(@times, mask_aa, ft_w) ));

% Advective term in physical space
wgradT = sum(bsxfun(@times,w_aa,gradb_aa),3);
% NB : in the stochastic case, w included both continuous and
% non-continuous components of the velocity

% Advective term in Fourier space
adv1 = - fft2(wgradT);clear wgradT

% % Remove aliasing
% adv1(ZM(1),:,:,:)=0;
% adv1(:,ZM(2),:,:)=0;
% %clear gradb

%% Laplacian diffusion term (from the stochastic material derivative)
%if ~isinf(model.sigma.k_c)
if model.sigma.sto
    adv2 = fct_heterogeneous_diff(model,ikx,iky,gradb,fft_b);
    %     if model.sigma.Smag.bool
    %         % Coefficient coef_Smag to target a specific diffusive scale
    %         adv2 = model.sigma.Smag.coef_Smag * adv2 ;
    %         % Taking into account the noise in the energy budget
    %         adv2 = model.advection.coef_diff * adv2;
    %     end
else
    adv2 = 0;
end
% if model.sigma.a0 > 0
%     if model.sigma.Smag.bool
% %         %%
% %         ikx = 1i*model.grid.k.kx;
% %         iky = 1i*model.grid.k.ky;
% %         gradb(:,:,1,:) = real(ifft2( ikx.*fft_b ));
% %         gradb(:,:,2,:) = real(ifft2( iky.*fft_b ));
%         adv2 = fct_heterogeneous_diff(model,ikx,iky,gradb,gradb_aa);
%         %         adv2 = fct_heterogeneous_diff(model,ikx,iky,fft_b,gradb,gradb_aa);
%         %         %%
%         %         %     adv2 = fct_heterogeneous_diff(model,ikx_aa,iky_aa,fft_b,gradb_aa);
%         %         %%
%         % Coefficient coef_Smag to target a specific diffusive scale
%         adv2 = model.sigma.Smag.coef_Smag * adv2 ;
%         % Taking into account the noise in the energy budget
%         adv2 = model.advection.coef_diff * adv2;
%     elseif model.sigma.hetero_modulation
%         adv2 = fct_heterogeneous_diff(model,ikx,iky,gradb,gradb_aa);
%
%     else
%         adv2 = - model.advection.coef_diff * k2 .* fft_b ;
%     end
% else
%     adv2 = 0;
% end

%% Hyperviscosity or turbulent dissipation
if model.advection.Lap_visco.bool | model.advection.HV.bool
    if model.advection.Smag.bool
        % Heterogeneous HV or diffusivity/viscosity coefficient
        if model.advection.Lap_visco.bool
            %             %%
            %             ikx = 1i*model.grid.k.kx;
            %             iky = 1i*model.grid.k.ky;
            %             gradb(:,:,1,:) = real(ifft2( ikx.*fft_b ));
            %             gradb(:,:,2,:) = real(ifft2( iky.*fft_b ));
            adv4 = fct_heterogeneous_diff(model,ikx,iky,gradb,fft_b);
            %             adv4 = fct_heterogeneous_diff(model,ikx,iky,fft_b,gradb);
            %             %%
            %             %adv4 = fct_heterogeneous_diff(model,ikx_aa,iky_aa,fft_b,gradb_aa);
            %             %%
            
            %             % Coefficient coef_Smag to target a specific diffusive scale
            %             adv4 = model.advection.Smag.coef_Smag * adv4 ;
            
            % Possibly add constant value
            % (but treated further below, without anti-aliasing)
            model.advection.HV.val = model.advection.HV.val * ...
                model.advection.Smag.weight_cst_dissip;
            %         coef_diff_aa = model.advection.HV.val * ...
            %             model.advection.HV.weight_cst_dissip ...
            %             + coef_diff_aa;
            
            
        elseif model.advection.HV.bool
            % Laplacian at the power p of buoyancy
            Lap_p_b_aa = bsxfun(@times, ...
                (-k2_aa) .^ (model.advection.HV.order/4) , ...
                fft_b );
            Lap_p_b_aa = real(ifft2(Lap_p_b_aa));
            
            %             % Heterogeneous HV coefficient
            %             coef_HV_aa = fct_coef_HV(model,fft_b,Lap_p_b_aa);
            
            % Possibly add constant value
            % (but treated further below, without anti-aliasing)
            model.advection.HV.val = model.advection.HV.val * ...
                model.advection.Smag.weight_cst_dissip;
            %         coef_HV_aa = model.advection.HV.val * ...
            %             model.advection.HV.weight_cst_dissip ...
            %             + coef_HV_aa;
            
            % HV coef * Laplacian at the power p of buoyancy
            adv4 = - model.advection.coef_diff .* Lap_p_b_aa;
            %             adv4 = - coef_HV_aa .* Lap_p_b_aa;
            
            adv4 = fft2(adv4);
            adv4 = bsxfun(@times, ...
                (-k2_aa) .^ (model.advection.HV.order/4) , ...
                adv4 );
            
        else
            error('Unknown deterministic subgrid tensor');
        end
        
    else
        adv4 = 0;
        
        %% Test
        %     adv4ref=adv4;
        %     %     % Norm of wave vector
        %     k2=model.grid.k.k2;
        %     %     k2=model.grid.k_HV.k2;
        %     %     warning('this k2 need to be smoother');
        %
        %
        %     if model.advection.Lap_visco.bool
        %         %% Test Laplacian
        %         % Visco/diff coef * Laplacian of buoyancy
        %         adv4 = bsxfun(@times, model.advection.HV.val, gradb);
        %         adv4 = fft2(adv4);
        %
        %         % Remove aliasing
        %         adv4(PX(1)+1,:,:,:)=0;
        %         adv4(:,PX(2)+1,:,:)=0;
        %
        %         % Divergence
        %         adv4 = 1i * kx .* adv4(:,:,1,:) + 1i * ky .* adv4(:,:,2,:);
        %
        %     elseif model.advection.HV.bool
        %         %% Test HV
        %         % Laplacian at the power p of buoyancy
        %         Lap_p_b = (-k2) .^ (model.advection.HV.order/4) .* fft_b;
        %         Lap_p_b = real(ifft2(Lap_p_b));
        %
        %         % HV coef * Laplacian at the power p of buoyancy
        %         adv4 = - model.advection.HV.val .* Lap_p_b;
        %         adv4 = fft2(adv4);
        %
        %         % Remove aliasing
        %         adv4(PX(1)+1,:,:,:)=0;
        %         adv4(:,PX(2)+1,:,:)=0;
        %
        %         % Laplacian at the power p
        %         adv4 = (-k2) .^ (model.advection.HV.order/4) .*  adv4;
        %     end
        %     %%
        %     err = abs(adv4-adv4ref)/mean(abs(adv4ref));
        %     mean(err)
        %
        %%
    end
    
    % Hyperviscosity term
    adv4 = adv4 - bsxfun( @times, ...
        model.advection.HV.val * k2 .^ (model.advection.HV.order/2) ,...
        fft_b);
else
    adv4 = 0;
end

%% Summing terms
d_fft_b_adv=adv1+adv2+adv4; clear adv1 adv2 adv4

% Remove aliasing
d_fft_b_adv(ZM(1),:,:,:)=0;
d_fft_b_adv(:,ZM(2),:,:)=0;

%% Forcing
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            d_fft_b_adv = bsxfun(@plus, d_fft_b_adv, ...
                +  model.advection.forcing.F );
        case 'Spring'
            d_fft_b_adv = d_fft_b_adv ...
                - model.advection.forcing.on_T * ...
                bsxfun(@plus, fft_b, - model.advection.forcing.F);
            %     d_fft_b_adv = d_fft_b_adv + model.advection.forcing.F;
    end
end

    function adv_hetero_diff = ...
            fct_heterogeneous_diff(model,ikx,iky,gradb,fft_b)
        %             fct_heterogeneous_diff(model,ikx,iky,gradb,gradb_aa)
        %         %            fct_heterogeneous_diff(model,ikx,iky,fft_b,gradb)
        % Compute heterogeneous Laplacian diffusion term
        %
        
        %         if (model.advection.Smag.bool & model.advection.Lap_visco.bool)
        % %         if (model.advection.Smag.bool & model.advection.Lap_visco.bool) ...
        % %                 | ( model.sigma.a0 > 0 & model.sigma.Smag.bool )
        %             % Heterogeneous dissipation coefficient
        %             coef_diff_aa = fct_coef_diff(model,nan,gradb_aa);
        %             % coef_diff_aa = fct_coef_diff(model,fft_b,gradb);
        %
        % %         elseif ( model.sigma.a0 > 0 & model.sigma.hetero_modulation)
        % %             coef_diff_aa = model.sigma.a0/2 * ...
        % %                 fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
        %         else
        %             coef_diff_aa = model.advection.coef_diff;
        %         end
        
        coef_diff_aa = model.advection.coef_diff;
        if model.advection.Smag.spatial_scheme
            b = real(ifft2(fft_b));
            gradb = gradient_mat_2_per(b,model.grid.dX);
            adv_hetero_diff = bsxfun(@times, coef_diff_aa, gradb);
            d_adv_dx = gradient_mat_2_per(adv_hetero_diff(:,:,1,:),model.grid.dX);
            d_adv_dy = gradient_mat_2_per(adv_hetero_diff(:,:,2,:),model.grid.dX);
            adv_hetero_diff = d_adv_dx(:,:,1,:) + d_adv_dy(:,:,2,:);
            adv_hetero_diff = fft2(adv_hetero_diff);
        else
            % Visco/diff coef * gradient of buoyancy
            adv_hetero_diff = bsxfun(@times, coef_diff_aa, gradb);
            adv_hetero_diff = fft2(adv_hetero_diff);
            
            % Divergence of ( Visco/diff coef * gradient of buoyancy )
            adv_hetero_diff = ikx .* adv_hetero_diff(:,:,1,:) ...
                + iky .* adv_hetero_diff(:,:,2,:);
        end
        %% Test
%         %adv_hetero_diff_ref = real(ifft2(adv_hetero_diff;
%         adv_hetero_diff_ref = adv_hetero_diff;
%         
%         b = real(ifft2(fft_b));
%         gradb = gradient_mat_2_per(b,model.grid.dX);
%         adv_hetero_diff = bsxfun(@times, coef_diff_aa, gradb);
%         d_adv_dx = gradient_mat_2_per(adv_hetero_diff(:,:,1,:),model.grid.dX);
%         d_adv_dy = gradient_mat_2_per(adv_hetero_diff(:,:,2,:),model.grid.dX);
%         adv_hetero_diff = d_adv_dx(:,:,1,:) + d_adv_dy(:,:,2,:);
%         adv_hetero_diff = fft2(adv_hetero_diff);
%         
%         adv_hetero_diff = real(ifft2(adv_hetero_diff));
%         adv_hetero_diff_ref = real(ifft2(adv_hetero_diff_ref));
%         err = abs(adv_hetero_diff - adv_hetero_diff_ref) / ...
%             mean(abs(adv_hetero_diff_ref(:)));
%         mean(err(:))
%             
%         figure(99);imagesc(((adv_hetero_diff')));axis equal;axis xy;colorbar;
%         figure(101);imagesc(((adv_hetero_diff_ref))');axis equal;axis xy;colorbar;
%         figure(104);imagesc(((err))');axis equal;axis xy;colorbar;
%         keyboard;
        %%
    end
end