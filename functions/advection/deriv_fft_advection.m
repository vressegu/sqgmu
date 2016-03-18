function d_fft_b_adv = deriv_fft_advection(model, fft_b, w)
% Compute the Fourier transform of the partial derivative along time 
% of the tracer according to the (stochastic or deterministic) transport
% equation with a hyperviscosity term
%

%% Grid of wave vectors
PX=model.grid.MX/2;
kx=model.grid.k.kx;
ky=model.grid.k.ky;
k2=model.grid.k.k2;

%% Advection term

% Gradient in Fourier space
adv1x = - 1i * kx .* fft_b;
adv1y = - 1i * ky .* fft_b;

% Gradient in physical space
gradb(:,:,1)=real(ifft2(adv1x));
gradb(:,:,2)=real(ifft2(adv1y));

% Advective term in physical space
wgradT=sum(bsxfun(@times,w,gradb),3);
% NB : in the stochastic case, w included both continuous and
% non-continuous components of the velocity

% Advective term in Fourier space
adv1=fft2(wgradT);clear wgradT

% Remove aliasing
adv1(PX(1)+1,:)=0;
adv1(:,PX(2)+1)=0;
clear gradb

%% Laplacian diffusion term (from the stochastic material derivative)
adv2 = - model.advection.coef_diff * k2 .* fft_b ;

%% Hyperviscosity

% Norm of wave vector
k2=model.grid.k_HV.k2;

% Hyperviscosity term
adv4 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;

%% Summing terms
d_fft_b_adv=adv1+adv2+adv4; clear adv1 adv2 adv4

