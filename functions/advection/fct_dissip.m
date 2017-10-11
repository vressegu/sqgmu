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
k2=model.grid.k_HV.k2;

%% Preprocessing

% Boolean filter

% Buoyancy
b = real(ifft2(fft_b));

% Gradient in Fourier space
adv1x = 1i * kx .* fft_b;
adv1y = 1i * ky .* fft_b;

% Gradient in physical space
gradb(:,:,1)=real(ifft2(adv1x));
gradb(:,:,2)=real(ifft2(adv1y));

% % Influence of the complex brownian variance
% sigma_on_sq_dt = sqrt(prod(model.grid.MX)) * sigma_on_sq_dt;

%% Advection term



if isinf(model.sigma.k_c) % Deterministic case
    noise_intake = zeros(size(b));
else % Stochastic case
    % Fourier transform of white noise
    dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 model.advection.N_ech]));
    % Multiplication by the Fourier transform of the kernel \tilde \sigma
    fft_sigma_dBt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
    clear dBt_C_on_sq_dt
    % Homogeneous velocity field
    sigma_dBt_dt = real(ifft2(fft_sigma_dBt));
    clear fft_sigma_dBt
    
    
    sigma_dBt_dtgradT2 = nan([model.grid.MX 1 model.advection.N_ech]);
    
    parfor sampl=1:model.advection.N_ech
        % Advective term in physical space
        sigma_dBt_dtgradT2(:,:,:,sampl)=sum(bsxfun(@times,sigma_dBt_dt(:,:,:,sampl),gradb),3);
        % NB : in the stochastic case, w included both continuous and
        % non-continuous components of the velocity
        
    end
        % Remove aliasing
%         sigma_dBt_dtgradT2(PX(1)+1,:,:,sampl)=0;
%         sigma_dBt_dtgradT21(:,PX(2)+1,:,sampl)=0;
    sigma_dBt_dtgradT2=fft2(sigma_dBt_dtgradT2);
    sigma_dBt_dtgradT2(PX(1)+1,:,:,:)=0;
    sigma_dBt_dtgradT2(:,PX(2)+1,:,:)=0;
    sigma_dBt_dtgradT2=real(ifft2(sigma_dBt_dtgradT2));
        
        
        %sigma_dBt_dtgradT2(:,:,:,sampl)=sigma_dBt_dtgradT2(:,:,:,sampl);
    %end
    
    noise_intake = model.advection.dt_adv/2 * mean(sigma_dBt_dtgradT2.^2,4);
    
    
end


%% Laplacian diffusion term (from the stochastic material derivative)
turb_dissip = model.advection.coef_diff * sum(gradb.^2,3);
% adv2 = - model.advection.coef_diff * k2 .* fft_b ;
% turb_dissip = - b .* real(ifft2(adv2));

%% Hyperviscosity

% adv4 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;
% dissip_HV = - b .* real(ifft2(adv4));
adv4 = k2 .^ (model.advection.HV.order/4) .* fft_b;
adv4 = real(ifft2(adv4));
dissip_HV = model.advection.HV.val * adv4.^2;

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


