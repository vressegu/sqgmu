function [coef_a_aa, coef_a] ...
    = fct_coef_estim_AbsDiff_heterogeneous(model,fft_w)
% Compute a heterogeneous multiplicative corrective coefficient to compute
% the heterogeneous aboslute diffusivity

if model.sigma.hetero_modulation
    % Velocity gradients
    fft_grad_v(:,:,:,:,1) = fct_grad(model,fft_w(:,:,1,:));
    fft_grad_v(:,:,:,:,2) = fct_grad(model,fft_w(:,:,2,:));
    fft_grad_v = permute(fft_grad_v,[1 2 5 4 3]);
    grad_v = real(ifft2(fft_grad_v)); clear fft_grad_v
    
    % Norm of velocity gradients
    n_grad_v = sqrt(sum(sum(grad_v.^2,5),3)); clear grad_v
end

% Velocity
v = real(ifft2(fft_w)); clear fft_w

% Square norm of the velocity
n_v2 = sum(v.^2,3); clear v

% mask_aa_LS = model.grid.k_aa.mask; %anti-aliasing mask
mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask

if model.sigma.hetero_modulation
    % Filtering the diffusivity coefficient at large scales
    n_v2 = fft2(n_v2 );
    n_v2 = real(ifft2( bsxfun(@times, n_v2, mask_aa_LS) ));
    
    % Filtering the diffusivity coefficient at large scales
    n_grad_v = fft2(n_grad_v);
    n_grad_v = real(ifft2( bsxfun(@times, n_grad_v, mask_aa_LS) ));
end

% Scaling of the absolute diffusivity
if model.sigma.hetero_modulation
    % V^2 / T ~ ||v||^2 / || grad(v) ||
    coef_a = n_v2 ./ (n_grad_v/sqrt(2)); %
elseif model.sigma.hetero_modulation_V2
    % V^2 ~ ||v||^2
    coef_a = n_v2; %
end
% abs_diff_w_estim_global1 = mean(abs_diff_estim_local(:))
% abs_diff_w_estim_global2 = mean(n_v2(:))/mean(n_grad_v(:)/sqrt(2))

% coef_a = 1/mean(abs_diff_estim_local(:)) * abs_diff_estim_local ;
% %coef_a = 1/model.sigma.a0 * abs_diff_estim_local ;

% Filtering the diffusivity coefficient at large scales
coef_a_aa = fft2(coef_a);
coef_a_aa = real(ifft2( bsxfun(@times, coef_a_aa, mask_aa_LS) ));

% Set negative values to zero
coef_a_aa = fct_smooth_pos_part(coef_a_aa);
% siz = size(coef_a_aa);
% coef_a_aa(coef_a_aa < 0)=0;
% coef_a_aa = reshape(coef_a_aa,siz);

% Normalisation
if nargout > 1
    coef_a = bsxfun(@times, 1/mean(mean(coef_a,2),1) , coef_a );
end
coef_a_aa = bsxfun(@times, 1/mean(mean(coef_a_aa,2),1) , coef_a_aa );

% % % warning('ft_w2 can be replaced by ft_w');
% % % warning('This treatement should be done realisation by realisation');
% %
% figure;imagesc(coef_a');colorbar;
% figure;imagesc(coef_a_aa');colorbar;
% % figure;imagesc(abs_diff_estim_local');
% % figure;imagesc(n_v2');
% % figure;imagesc(n_grad_v');
% figure;imagesc(coef_a');
