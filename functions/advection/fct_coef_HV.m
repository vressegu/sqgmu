function [coef_HV_aa,coef_HV] = fct_coef_HV(model,fft_b,Lap_p_b_aa)
% Compute a heterogeneous HV coefficient which target a specific
% dissipation scale

%% Filter
mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask

if nargin < 3
    %% Grid of wave vectors
    k2_aa = model.grid.k_aa.k2;
    
    %% Laplacian at the power p of buoyancy  anti-aliased
    Lap_p_b_aa = (-k2_aa) .^ (model.advection.HV.order/4) .* fft_b;
    Lap_p_b_aa = real(ifft2(Lap_p_b_aa));
end

coef_HV = sqrt(abs( Lap_p_b_aa ));

% Coefficient coef_Smag to target a specific diffusive scale
coef_HV = model.advection.Smag.coef_Smag * coef_HV ;

% Filtering the HV coefficient at large scales
coef_HV = fft2(coef_HV );
coef_HV_aa = real(ifft2( bsxfun(@times, coef_HV, mask_aa_LS) ));

