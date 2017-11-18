function nrj = fct_mean_noise_energy(model,fft_b)
% Compute the spatial mean of the energy intake of the noise
%

% %% Grid of wave vectors
% PX=model.grid.MX/2;
% kx=model.grid.k.kx;
% ky=model.grid.k.ky;
% k2=model.grid.k.k2;
% 
% % %% Advection term
% % 
% % % Gradient in Fourier space
% % adv1x = 1i * kx .* fft_b;
% % adv1y = 1i * ky .* fft_b;
% % 
% % % Gradient in physical space
% % gradb(:,:,1)=ifft2(adv1x);
% % gradb(:,:,2)=ifft2(adv1y);

%% Energy intake

% Kernel
kern = fct_integral_noise_energy(model,model.grid.k.kx,model.grid.k.ky);

% Energy intake of the noise by wave-vector
nrj = abs(fft_b).^2 .* kern;

% Integration over the wave-vectors
% nrj = nrj(:);
% len = length(nrj);
% nrj = sum(nrj);
nrj = sum(sum(nrj,1),2);

nrj = nrj / prod(model.grid.MX); % Discrete FFT ?
% nrj = nrj / len; % Discrete FFT ?

