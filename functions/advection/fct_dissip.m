function [turb_dissip, noise_intake, dissip_HV, intake_forcing] = ...
    fct_dissip(model,fft_b,sigma_on_sq_dt)
% Compute the deterministic and stochastic dissipation terms locally in
% space.
%

warning('This function is only valid for homogeneous dissipation');



%% Grid of wave vectors
PX=model.grid.MX/2;
kx=model.grid.k.kx;
ky=model.grid.k.ky;
k2=model.grid.k.k2;
kx1 = permute(kx,[ 3 4 1 2]);
ky1 = permute(ky,[ 3 4 1 2]);

% Boolean filter

% Buoyancy
b = real(ifft2(fft_b));

%% Advection term

% warning('pk un signe - ici ?');
% % Gradient in Fourier space
% adv1x = - 1i * kx .* fft_b;
% adv1y = - 1i * ky .* fft_b;
% 
% % Influence of the complex brownian variance
% sigma_on_sq_dt = sqrt(prod(MX)) * sigma_on_sq_dt;
% 
% % Gradient in physical space
% gradb(:,:,1)=real(ifft2(adv1x));
% gradb(:,:,2)=real(ifft2(adv1y));
% 
% % % Advective term in physical space
% % wgradT=sum(bsxfun(@times,w,gradb),3);
% % % NB : in the stochastic case, w included both continuous and
% % % non-continuous components of the velocity
% % 
% % % Advective term in Fourier space
% % adv1=fft2(wgradT);clear wgradT
% % 
% % % Remove aliasing
% % adv1(PX(1)+1,:)=0;
% % adv1(:,PX(2)+1)=0;
% % clear gradb

noise_intake = zeros(size(b));

%% Laplacian diffusion term (from the stochastic material derivative)
adv2 = - model.advection.coef_diff * k2 .* fft_b ;
turb_dissip = - b .* real(ifft2(adv2));

%% Hyperviscosity

% Norm of wave vector
k2=model.grid.k_HV.k2;

% Hyperviscosity term
adv4 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;
dissip_HV = - b .* real(ifft2(adv4));

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


