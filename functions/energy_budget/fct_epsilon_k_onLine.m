%function [PI_loc_smoothKX_pos_n,PI_loc_smoothKX,PI_loc_smoothK_neg,PI_loc_smoothK]...
function [coef_AbsDiff, PI_loc_smoothKX_pos] ...
    = fct_epsilon_k_onLine(model,fft_b,fft_w,fft_gradb_aa)
%     = fct_epsilon_k_onLine(model,fft_b,w_aa,fft_gradb_aa)
% Compute the spatially-heterogeneous energy flux accross scale averaged
% over some wave numbers.
% The result is normalised in order to have a spatial mean equal to 1.
%

bool_plot = false;

if model.advection.N_ech > 1
    %    warning('This function is not optimized for ensemble forecasts');
end

%% Grid of wave vectors
ZM = model.grid.k.ZM; %index of modes ot be zero'ed out
mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask
mask_aa = model.grid.k_aa.mask; %anti-aliasing mask

M_kappa=min(model.grid.MX);
P_kappa= M_kappa/2;
d_kappa = 2*pi/sqrt(prod(model.grid.MX.* model.grid.dX));
kappa= d_kappa * ( 0:(P_kappa-1) ) ;

k = model.grid.k.k;

%% Filtering of intial buoyancy
if bool_plot
    figure(55);imagesc(real(ifft2(fft_b))');axis xy;axis equal
end
if model.sigma.hetero_energy_flux_prefilter
    fft_b = bsxfun(@times, ...
        model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
        fft_b);
end
if bool_plot
    figure(56);imagesc(real(ifft2(fft_b))');axis xy;axis equal
end

%% Gradient of b, anti-aliased
% if nargin < 3
fft_w = SQG_large_UQ(model, fft_b);
% end
% if nargin < 4
ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
iky_aa = model.grid.k_aa.iky;
% in Fourier space, de-aliased, then in physical space.
fft_gradb_aa(:,:,1,:) = bsxfun(@times,ikx_aa,fft_b);
fft_gradb_aa(:,:,2,:) = bsxfun(@times,iky_aa,fft_b);
gradb_aa =  real(ifft2(fft_gradb_aa ));
% end

%% dealisasing of the velocity: FT, apply mask, iFT
% ft_w = fft2( w );
w_aa = real(ifft2( bsxfun(@times, mask_aa, fft_w) ));

%% Advection term

% Advective term in physical space
wgradT = sum(bsxfun(@times,w_aa,gradb_aa),3);
%         wgradT = sum(bsxfun(@times,w_aa(:,:,:,sampl),gradb_aa_SS),3);
fft_wgradT = fft2(wgradT);

%%
PI_loc_smoothK = zeros([model.grid.MX 1 model.advection.N_ech]);
for sampl=1:model.advection.N_ech
    kappa_trunc = ( model.sigma.km_LS(sampl) + d_kappa) : d_kappa : ...
        (model.sigma.kappaLSforEspi_on_kappamin * model.sigma.kappamin_on_kappamax * kappa(end) - d_kappa);
    %     kappa_trunc = ( model.sigma.km_LS(sampl) + d_kappa) : d_kappa : ...
    %         (model.sigma.kappamin_on_kappamax * kappa(end) - d_kappa);
    for  kappa_local = kappa_trunc
        %% Separating large scales and small-scales
        
        % Filter at scale kappa
        alpha = 1.;
        order = 19.;
        steep_filter = exp(-alpha*( 1/(eps+kappa_local) ...
            .* k ).^order );
        steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
        steep_filter(:,ZM(2)) = 0.;
        
        % Buoyancy fitlered at scale kappa
        b_LS = real(ifft2( steep_filter .* fft_b(:,:,:,sampl) ));
        
        % Buoyancy gradient high-pass filtered at scales kappa
        %         fft_gradb_aa_SS = fft_gradb_aa(:,:,:,sampl) ;
        % %         warning('High-pass filter removed')
        %         %         fft_gradb_aa_SS = bsxfun( @times, (1 - steep_filter) , ...
        %         %             fft_gradb_aa(:,:,:,sampl) ) ;
        %         fft_gradb_aa_SS(ZM(1),:,:,:) = 0.; %de-alias the single high freq
        %         fft_gradb_aa_SS(:,ZM(2),:,:) = 0.;
        %         gradb_aa_SS =  real(ifft2(fft_gradb_aa_SS ));
        
        %% Advection term
        
        %         % Advective term in physical space
        %         wgradT = sum(bsxfun(@times,w_aa(:,:,:,sampl),gradb_aa),3);
        % %         wgradT = sum(bsxfun(@times,w_aa(:,:,:,sampl),gradb_aa_SS),3);
        %         wgradT = real(ifft2( steep_filter .* fft2(wgradT) ));
        
        wgradT = real(ifft2( steep_filter .* fft_wgradT(:,:,:,sampl) ));
        
        % Spatially local flux at scale kappa
        PI_loc_smoothK(:,:,:,sampl) = PI_loc_smoothK(:,:,:,sampl) ...
            + bsxfun( @times, b_LS , wgradT);
        
    end
    % Averaging over spatial frequencies
    PI_loc_smoothK(:,:,:,sampl) = PI_loc_smoothK(:,:,:,sampl) ...
        / length(kappa_trunc);
end
% Desaliasing
PI_loc_smoothK = fft2(PI_loc_smoothK);

% % if model.sigma.hetero_energy_flux_postfilter
%     PI_loc_smoothK = bsxfun(@times, ...
%         model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
%         PI_loc_smoothK);
% % end

PI_loc_smoothK(ZM(1),:,:,:) = 0.; % de-alias the single high freq
PI_loc_smoothK(:,ZM(2),:,:) = 0.; % de-alias the single high freq
PI_loc_smoothK = real(ifft2(PI_loc_smoothK));

%%

% if model.sigma.hetero_energy_flux_v2
%     if bool_plot
%         figure(42);imagesc(PI_loc_smoothK');axis xy;axis equal
%     end
%     PI_loc_smoothK = fct_smooth_pos_part(PI_loc_smoothK);
%     PI_loc_smoothK = PI_loc_smoothK .^ (1/3) ;
% end

%%

if bool_plot
    figure(48);imagesc(PI_loc_smoothK');axis xy;axis equal
end

if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
        &&    model.sigma.hetero_energy_flux_postfilter
    mask_post_filter = model.grid.k_aa_sigma_hetero_energy_flux_post.mask;
    % mask_post_filter = mask_aa_LS;
    PI_loc_smoothK = real(ifft2( bsxfun(@times, fft2(PI_loc_smoothK), ...
        mask_post_filter) ));
    % PI_loc_smoothKX = real(ifft2( bsxfun(@times, PI_loc_smoothKX, ...
    %     mask_aa_LS) ));
end

if bool_plot
    figure(49);imagesc(((PI_loc_smoothK))');axis xy;axis equal
end

%%
% if ~ model.sigma.hetero_energy_flux_v2
%%
% figure(88);imagesc(PI_loc_smoothK');axis xy;axis equal
% PI_loc_smoothK = real(ifft2( bsxfun(@times, VLS_steep_filter, ...
%     fft2(PI_loc_smoothK) ) ));
% figure(89);imagesc(PI_loc_smoothK');axis xy;axis equal

%% Removing negative value of the energy transfer, but smoothing the
% negative values in order to take into account their influence
% in the global budget

% Filter at very large scales

if model.sigma.hetero_energy_flux_v2
    kappa_VLS = model.sigma.kappa_VLS_on_kappa_LS ...
        * model.sigma.kappamin_on_kappamax * kappa(end);
else
    kappa_VLS = model.sigma.kappa_VLS_on_kappa_LS ...
        * model.advection.Smag.dealias_ratio_mask_LS * kappa(end);
end
% kappa_VLS = 1/4 * model.advection.Smag.dealias_ratio_mask_LS * kappa(end);
alpha = 1.;
order = 19.;
VLS_steep_filter = exp(-alpha*( 1/kappa_VLS .* k ).^order );
VLS_steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
VLS_steep_filter(:,ZM(2)) = 0.;

%%

% Negative part
PI_loc_smoothK_neg = - fct_smooth_pos_part( - PI_loc_smoothK );

% Smoothing of the negative part at very large scales
PI_loc_smoothK_neg = real(ifft2( bsxfun(@times, VLS_steep_filter, ...
    fft2(PI_loc_smoothK_neg) ) ));

% Positive part
PI_loc_smoothK_pos = fct_smooth_pos_part(PI_loc_smoothK);

% Gather postive and filtered negative part
PI_loc_smoothK_pre_pos = PI_loc_smoothK_pos + PI_loc_smoothK_neg;


if bool_plot
    figure(28);imagesc(((PI_loc_smoothK_pre_pos))');axis xy;axis equal
    %     close all;
end


%% Filtering the coefficient at large scales
PI_loc_smoothKX = fft2(PI_loc_smoothK_pre_pos);


%%
%
% figure(48);imagesc(real(ifft2(PI_loc_smoothKX))');axis xy;axis equal
%
% mask_post_filter = model.grid.k_aa_sigma_hetero_energy_flux.mask;
% mask_post_filter = mask_aa_LS;
% PI_loc_smoothKX = real(ifft2( bsxfun(@times, PI_loc_smoothKX, ...
%     mask_post_filter) ));
% % PI_loc_smoothKX = real(ifft2( bsxfun(@times, PI_loc_smoothKX, ...
% %     mask_aa_LS) ));
%
% figure(49);imagesc(((PI_loc_smoothKX))');axis xy;axis equal
%%

PI_loc_smoothKX = real(ifft2( bsxfun(@times, PI_loc_smoothKX, ...
    mask_aa_LS) ));

% Forced positiveness
%     PI_loc_smoothKX_pos = fct_smooth_pos_part(PI_loc_smoothKX);
%
%     % %% Coefficient for the Absolute diffusivity by scale
%     % coef_AbsDiff = PI_loc_smoothKX_pos .^ (1/3) ;
%     %
%     % %% Normalisation
%     % coef_AbsDiff = bsxfun(@times, 1./mean(mean(coef_AbsDiff,2),1), ...
%     %     coef_AbsDiff );


% else
%     PI_loc_smoothKX = PI_loc_smoothK;
% end
% Forced positiveness
PI_loc_smoothKX_pos = fct_smooth_pos_part(PI_loc_smoothKX);

if bool_plot
    figure(50);imagesc(((PI_loc_smoothKX_pos))');axis xy;axis equal
    %     close all;
end

%% Normalisation
if model.sigma.hetero_energy_flux_averaging_after
    PI_loc_smoothKX_pos = PI_loc_smoothKX_pos .^ (1/3) ;
end
coef_AbsDiff = bsxfun(@times, 1./mean(mean(PI_loc_smoothKX_pos,2),1), ...
    PI_loc_smoothKX_pos );

%% Coefficient for the Absolute diffusivity by scale
if ~ model.sigma.hetero_energy_flux_averaging_after
    coef_AbsDiff = coef_AbsDiff .^ (1/3) ;
end


if bool_plot
    figure(66);imagesc(((coef_AbsDiff))');axis xy;axis equal
    close all;
end

