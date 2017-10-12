function fft_grad_f = fct_grad(model,fft_f)
% Compute the gradient of the function f in Fourier space
%

%% Grid of wave vectors
kx=model.grid.k.kx;
ky=model.grid.k.ky;

%% Gradient in Fourier space
fft_grad_f(:,:,1,:) = 1i * kx .* fft_f;
fft_grad_f(:,:,2,:) = 1i * ky .* fft_f;

