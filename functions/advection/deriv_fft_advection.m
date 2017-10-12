function d_fft_b_adv = deriv_fft_advection(model, fft_b, w)
% Compute the Fourier transform of the partial derivative along time
% of the tracer according to the (stochastic or deterministic) transport
% equation with a hyperviscosity term
%

%% Grid of wave vectors
PX=model.grid.MX/2;
% Stable k
kx=model.grid.k.kx;
ky=model.grid.k.ky;
% Stable k2
k2=model.grid.k.k2;

%% Advection term

% Gradient in Fourier space
adv1x = 1i * kx .* fft_b;
adv1y = 1i * ky .* fft_b;

% Gradient in physical space
gradb(:,:,1,:)=real(ifft2(adv1x));
gradb(:,:,2,:)=real(ifft2(adv1y));

% Advective term in physical space
wgradT=sum(bsxfun(@times,w,gradb),3);
% NB : in the stochastic case, w included both continuous and
% non-continuous components of the velocity

% Advective term in Fourier space
adv1 = - fft2(wgradT);clear wgradT

% Remove aliasing
adv1(PX(1)+1,:,:,:)=0;
adv1(:,PX(2)+1,:,:)=0;
%clear gradb

%% Laplacian diffusion term (from the stochastic material derivative)
adv2 = - model.advection.coef_diff * k2 .* fft_b ;

%% Hyperviscosity

if model.advection.Smag.bool
    % Stable k2
    k2=model.grid.k.k2;
    
    % Heterogeneous HV or diffusivity/viscosity coefficient
    if model.advection.Lap_visco.bool
        % Root square of the buoyancy gradient norm
        model.advection.HV.val = sum(gradb.^2,3) .^(1/4) ;
        
        % Coefficient coef_Smag to target a specific diffusive scale
        model.advection.HV.val = model.advection.Smag.coef_Smag * ...
            model.advection.HV.val ;
        
        % Visco/diff coef * Laplacian of buoyancy
        adv4 = bsxfun(@times, model.advection.HV.val, gradb);
        adv4 = fft2(adv4);
        
        % Remove aliasing
        adv4(PX(1)+1,:,:,:)=0;
        adv4(:,PX(2)+1,:,:)=0;
        
        % Divergence
        adv4 = 1i * kx .* adv4(:,:,1,:) + 1i * ky .* adv4(:,:,2,:);
        
    elseif model.advection.HV.bool
        % Laplacian at the power p of buoyancy
        Lap_p_b = (-k2) .^ (model.advection.HV.order/4) .* fft_b;
        Lap_p_b = real(ifft2(Lap_p_b));
%         warning('this k2 need to be smoother ?');
        model.advection.HV.val = sqrt(abs( Lap_p_b ));
        
        % Coefficient coef_Smag to target a specific diffusive scale
        model.advection.HV.val = model.advection.Smag.coef_Smag * ...
            model.advection.HV.val ;
        
        % HV coef * Laplacian at the power p of buoyancy
        adv4 = - model.advection.HV.val .* Lap_p_b;
        adv4 = fft2(adv4);
        
        % Remove aliasing
        adv4(PX(1)+1,:,:,:)=0;
        adv4(:,PX(2)+1,:,:)=0;
        
        % Laplacian at the power p
        adv4 = (-k2) .^ (model.advection.HV.order/4) .*  adv4;
    else
        error('Unknown deterministic subgrid tensor');
    end
    
else
    % Norm of (unstable) wave vector
    k2=model.grid.k_HV.k2;
    
    % Hyperviscosity term
    adv4 = - model.advection.HV.val * k2 .^ ...
        (model.advection.HV.order/2) .* fft_b;
    
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

%% Summing terms
d_fft_b_adv=adv1+adv2+adv4; clear adv1 adv2 adv4

%% Forcing
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            d_fft_b_adv = d_fft_b_adv ...
                +  model.advection.forcing.F;
        case 'Spring'
            d_fft_b_adv = d_fft_b_adv ...
                - model.advection.forcing.on_T * ...
                ( fft_b - model.advection.forcing.F);
            %     d_fft_b_adv = d_fft_b_adv + model.advection.forcing.F;
    end
end
