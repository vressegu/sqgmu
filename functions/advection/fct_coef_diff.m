function [coef_diff_aa,coef_diff] = fct_coef_diff(model,fft_b,gradb_aa)
% Compute a heterogeneous dissipation coefficient which target a specific
% dissipation scale

%% Filter
mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask

if nargin < 3
    %% Grid of wave vectors
    ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
    iky_aa = model.grid.k_aa.iky;
    
    %% Gradient of b, anti-aliased
    % in Fourier space, de-aliased, then in physical space.
    gradb_aa(:,:,1,:) = real(ifft2( ikx_aa.*fft_b ));
    gradb_aa(:,:,2,:) = real(ifft2( iky_aa.*fft_b ));
    
end

% Root square of the buoyancy gradient norm
coef_diff = sum(gradb_aa.^2,3) .^(1/4) ;

% % Coefficient coef_Smag to target a specific diffusive scale
% coef_diff = model.advection.Smag.coef_Smag * coef_diff ;

% Filtering the diffusivity coefficient at large scales
coef_diff = fft2(coef_diff );
coef_diff_aa = real(ifft2( bsxfun(@times, coef_diff, mask_aa_LS) ));

