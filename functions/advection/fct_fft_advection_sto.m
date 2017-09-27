 function [fft_buoy_part, model] = fct_fft_advection_sto(model,  fft_buoy_part)
% Advection of buoyancy using SQG or SQG_MU model
%
% Modified by P. DERIAN 2016-08-18
%   - changed stochastic/deterministic tests;
%   - added support for several noise models (switches);
%   - improved verbose mode;
%   - improved output path generation for cross-platform compatibility;
%   - merged diffusion coefficient with noise models pre-computations;
%   - removed model.advection.coef_diff for homogeneization purpose;
%   - renamed "sigma_dBt_dt" as "sigma_dBt_over_dt" for clarity.
% Modified by P. DERIAN 2016-10-11
%   - made plotting optional, enabled custom output root directory.
% Modified by P. DERIAN 2016-10-[24--28]
%   - function flow changed to boost parallelization;
%   - moved sigma_on_sq_dt to model.sigma to comply with parfor;
%   - adjusting dt so that a day (86400 s) is a multiple of dt; 
%   - created a (nested) output() function for graphics and data;
% Modified by P. DERIAN 2017-05-04
%   - ensemble buoyancy initialization is now being done before entering
%       this function.
% TODO:
%   - output frequency as a model parameter;
%   - allow user to manually set the dt_adv;
%   - deal with restart (ouch...?);
%   - ...

%% Folder to save plots and files
% Generate the case name (initial condition + forcing if any)
config_name = sprintf('%s_%d', model.init.label_data, model.grid.MX(1));
if ~strcmp(model.advection.forcing.type_forcing, 'None') % Append name of forcing
    config_name = [config_name, '_', model.advection.forcing.type_forcing];
end
% Then the path, depending whether sotchastic or not
if model.is_stochastic % Stochastic case
    model.output.folder_simu = fullfile(model.output.folder_root, 'SQG_MU', ...
                                        config_name, model.sigma.type_noise);
else % Deterministic case  
    model.output.folder_simu = fullfile(model.output.folder_root, 'usual_SQG', ...
                                        config_name);
end
% Create the folders
fct_create_output_folders(model);
% Colormap
tmp = load('BuYlRd.mat');
model.output.colormap = tmp.BuYlRd; clear tmp;
% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.output.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid
% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% Grid in Fourier space
model = init_grid_k(model);
% Dealias the single high-freq of initial fft_b
fft_buoy_part(model.grid.k.ZM(1),:,:) = 0.;
fft_buoy_part(:,model.grid.k.ZM(2),:) = 0.;

%% Initialisation of the spatial velocity fields
N_ech = model.advection.N_ech; % ensemble size
% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_buoy_part);
w = real(ifft2(fft_w));
% [WIP] Forcing (stored in model)
model.advection.forcing.F = forcing(model); 
% This temporary field is for setting hyperviscosity and CFL
% [TODO] make it better?
tmp_w = mean(w,4) + model.advection.forcing.F;

%% Diffusion coeff
% [TODO] move within spectral... bound for CFL "bound1"?
model.sigma.a0 = 2.*model.physical_constant.f0/model.sigma.k_c^2; 

%% Hyperviscosity
% Root mean square Lyapunov
[dxUx,dyUx] = gradient(tmp_w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(tmp_w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s = (dxUx+dyUy)/2;
d = (dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s d dxUx dxUy dyUx dyUy lambda
% Order of the hyperviscosity
model.advection.HV.order = 8;
% Hyperviscosity coefficient
model.advection.HV.val = ...
    40 * model.advection.lambda_RMS * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;

%% Choice of time step : CFL
% CFL of the diffusion (or CFL of the white noise advection)
% [TODO] what if there is no TLU diffusion? (cf WavHypervis noise model)
dX2 = (model.grid.dX /pi).^2;
bound1 = 2./model.sigma.a0*prod(dX2)/sum(dX2);
% CFL of the (large-scale) advection
dX = permute(model.grid.dX,[1 3 2]);
bound2 = sum(bsxfun(@times,abs(tmp_w),pi./dX),3);
bound2 = max(bound2(:));
bound2 = 1./bound2/4.;
% CFL of the hyperviscosity
bound3 = 1./model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
         (model.advection.HV.order/2);
clear dX dX2
% Minimum of the CFL
dt = min([bound1 bound2 bound3]);
clear bound1 bound2 bound3
if model.is_stochastic
    dt = dt/2.;
    % Further constraint on dt due to the use of a (simple) Euler scheme 
    % for the SPDE
end
% [WIP] adjust dt so that a day (86400 s) is a multiple of dt
dt = floor(dt);
while mod(86400, dt)>0 && dt>2 % brute force...
    dt = dt - 1;
end
model.advection.dt_adv = double(dt);
clear dt tmp_w;
% Number of time step
N_t = ceil(model.advection.advection_duration/model.advection.dt_adv);

%% Pre-computations for noise models
switch model.sigma.type_noise
    case 'Spectrum'
        % Fourier transform of the kernel \tilde sigma up to a multiplicative
        % constant
        [sigma_on_sq_dt, ~, missed_var_small_scale_spectrum ] ...
            = fct_sigma_spectrum(model,fft_w);
        % sigma_on_sq_dt will be used to simulate sigma d B_t
        % missed_var_small_scale_spectrum will be used to set the mulitplicative
        % constant
        % Muliplicative constant of the kernel \tilde sigma
        model.sigma.a0_on_dt = model.sigma.a0/model.advection.dt_adv;
        model.sigma.sigma_on_sq_dt = sqrt(2.*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
            * sigma_on_sq_dt;
        % the factor d=2 is for the dimension d of the space R^d
        % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
        clear missed_var_small_scale_spectrum
    case 'SVDfull'
        model.sigma.noise_generator = SVDnoise3D(model.grid.MX, model.sigma.P, ...
                                                 model.sigma.boundary_condition);
        % Scaling factor for the basis
        % to compute a_l from large scale variance E[(u_L - E[u_L])^2]
        % Note: * sqrt(P^(-2/3)) is an amplitude scaling for the change of
        %           scale (observation=>noise);
        %       * the time scaling is computed on the fly later.
        model.sigma.scaling = model.sigma.noise_generator.default_amplitude_scaling();
        % Periodic gradient operators (in the physical space)
        % these are centered O(2) finite differences with periodic BC,
        % used to compute div(a).
        [D1, D2] = gradient_periodic(model.grid.MX(1), model.grid.MX(2), ...
                                     model.grid.dX(1), model.grid.dX(2));
        model.sigma.D1 = D1;
        model.sigma.D2 = D2;
        % Clean
        clear D1 D2
    case 'Hypervis' %[WIP]
        % The Scrambler class is used to generate pseudo-observations from local
        % permutations of a field (same as what happens in SVDnoise3d)
        model.sigma.scrambler = Scrambler(model.grid.MX, model.sigma.P, ...
                                          model.sigma.boundary_condition);
        % This is the stochastic advection operator (sigma_dBt . grad(b))
        % NB: b is de-aliased in the operator (model.grid.k_aa.mask)
        model.sigma.fft_StoAdv = sqrt(2.*model.sigma.scaling*model.advection.HV.val) * ...
                                (-model.grid.k.k2).^(model.advection.HV.order/4) .* ...
                                model.grid.k_aa.mask; 
    case 'WavHypervis' %[WIP]
        % The coarse wavelet scale
        model.sigma.Jmin = max([0, nextpow2(min(model.grid.MX))-model.sigma.N_lvl]);
        % The wavelet filter
        switch lower(model.sigma.wavelet_type)
            case 'daubechies'
                model.sigma.qmf = MakeONFilter('Daubechies', 2*model.sigma.wavelet_vm*2);
            case 'haar'
                model.sigma.qmf = MakeONFilter('Haar');
            case 'symmlet'
                model.sigma.qmf = MakeONFilter('Symmlet', model.sigma.wavelet_vm);
            otherwise
                error('SQGMU:fct_fft_advection_sto:ValueError', ...
                      'Wavelet type "%s" is currently not supported by the WavHypervis noise model.', ...
                      model.sigma.wavelet_type)
        end
        % This is the stochastic advection operator (sigma_dBt . grad(b))
        % NB: - sqrt(dt_adv) coeff for the BM variance.
        %     - no minus applied to model.advection.HV.val as it is
        %     positive in this code
        %     - HV.order/4 as the HV operator uses order/2, see deriv_fft_advection(). 
        model.sigma.fft_StoAdv = sqrt(2.*model.sigma.scaling*model.advection.HV.val*model.advection.dt_adv) * ...
                                (-model.grid.k.k2).^(model.advection.HV.order/4);
    case 'None'
        % nothing to do here.
    otherwise
        error('SQGMU:fft_advection_sto:ValueError', '"%s" is not a valid noise model', ...
              model.sigma.type_noise);
end

%% Allocations for output variables
% Note: w was created before, see "Initialisation of the spatial fields"
w = zeros([model.grid.MX 2 N_ech]); % time-correlated velocity
sigma_dBt_over_dt = zeros([model.grid.MX 2 N_ech]); % time-uncorrelated velocity
switch model.sigma.type_noise
    case {'SVDfull'}
        A = zeros([model.grid.MX 3 N_ech]); % Diffusion tensor data
    case {'Hypervis', 'WavHypervis'}
        sigma_dBt_gradb = zeros([model.grid.MX N_ech]); % Stochastic advection
end
% output frequency in time-step
% [TODO] set as a model parameter
output_frequency = ceil(86400/model.advection.dt_adv);
model.output.frequency = output_frequency;

%% [WIP] Time probes
% this is an experimental feature that probes the given variables every
% time-step.
has_probes = false;
probes = cell(N_ech,1);
if has_probes
    iProbe = [16, 32, 48, 64, 80, 96, 112, 128]; % probe indices
    jProbe = [16, 32, 48, 64, 80, 96, 112, 128];
    maxSize = N_t;
    for n=1:N_ech
        probes{n} = ProbeGrid(iProbe, jProbe, maxSize, {'u', 'v', 'b'});
    end
end

%% Display information
% Print some information
if model.verbose
    if model.is_stochastic
        fprintf(1, '\nStochastic simulation: %s ("%s")\n', ...
                model.sigma.long_name, model.sigma.type_noise);
    else
        fprintf(1, 'Deterministic simulation\n');
    end
    fprintf(1, 'Initial condition: %s (%s)\n', model.init.type_data, model.init.label_data);
    fprintf(1, 'Velocity forcing: %s\n', model.advection.forcing.type_forcing);
    fprintf(1, 'Output folders: %s\n', model.output.folder_simu);
    fprintf(1, 'Output frequency: %d steps (%.3f days)\n', model.output.frequency, ...
            model.output.frequency*model.advection.dt_adv/(3600*24.));
    if ~model.output.plot_results
        fprintf(1, '\t[!] Plotting disabled (-nojvm compatible mode)\n');
    end
    fprintf(1, '1/k_c is equal to %.2f m\n', 1./model.sigma.k_c);
    fprintf(1, 'Time step: %.2f seconds\n', model.advection.dt_adv);
    fprintf(1, 'Total time of advection: %f days\n', N_t*model.advection.dt_adv/3600./24.);
    fprintf(1, 'Grid size: %dx%d (de-aliasing="%s")\n', model.grid.MX(1), model.grid.MX(2), model.grid.dealias_method);
    fprintf(1, 'Ensemble size: %d realization(s)\n', N_ech);
    fprintf(1, 'Random generator seed: %d\n', model.seed);
    if has_probes
       fprintf(1, '\n[WIP] Warning! Time probes enabled...\n');
    end;
    fprintf(1, '\nBeginning simulation...\n\n');
end 

%% Parallel Business
try
    pool = parpool(); % create a pool of workers for parallel computation
                      % (as many workers as cores available)
catch
    fprintf(1, '\n[!] unable to start parallel pool. Is a pool running already?.\n');
    fprintf(1, '    Consider compiling the SQGMU software (see /compile) if Distributed Computing licenses are unavailable.\n\n');
end

%% Loop over time
t = 1; % [TODO] deal with restart
while t<=N_t   
    % write the initial state
    % [TODO] apply SQG_large_UQ() to save initial velocities?
    if t==1
        output()
    end
    
    % set the current time interval
    begin_t = t;
    end_t = min((t+model.output.frequency-1), N_t);
    fprintf(1 ,'%s - beginning iterations: [%d, %d]\n', datestr(datetime), begin_t, end_t);
    
    %% for each sample (s_)
    parfor sampl=1:N_ech
        % Note: define these dummy variables below just to suppress some PARFOR warning
        s_sigma_dBt_over_dt = 0.;
        s_sigma_dBt_gradb = 0.; %[WIP] only used by 'Hypervis' noise model 
        s_fft_sigma_dBt_gradb = 0.; %[WIP] only used by 'Hypervis' noise model 
        s_A = 0.;
        s_probe = 0;
        % retrieve the current buoyancy
        s_fft_b = fft_buoy_part(:,:,sampl);
        % [WIP] retrieve sample's probe
        if has_probes
            s_probe = probes{sampl};
        end;
        % iterate over the time interval
        for tt=begin_t:end_t
            
            %% Time-correlated velocity
            s_fft_w = SQG_large_UQ(model, s_fft_b);
            s_w = real(ifft2( s_fft_w )) + model.advection.forcing.F; %include optional forcing
            
            %% [WIP] update probe
            if has_probes
                s_probe.update_field('u', s_w(:,:,1));
                s_probe.update_field('v', s_w(:,:,2));
                %Note: it's annoying that it costs another FT for the
                %buoyancy
                s_probe.update_field('b', real(fft2(s_fft_b)));
            end
                
            %% Time-uncorrelated velocity 
            if model.is_stochastic
                switch model.sigma.type_noise
                    % Spectral noise (isotropic and homogeneous in space)
                    case 'Spectrum'
                        % Fourier transform of white noise
                        tmp_dBt_C_on_sq_dt = fft2( randn( model.grid.MX ));
                        % Multiplication by the Fourier transform of the kernel \tilde \sigma
                        tmp_fft_sigma_dBt = bsxfun(@times, ...
                                                   model.sigma.sigma_on_sq_dt, ...
                                                   tmp_dBt_C_on_sq_dt);
                        % Homogeneous velocity field
                        s_sigma_dBt_over_dt = real(ifft2(tmp_fft_sigma_dBt));
                        % Stochastic diffusion (scalar here)
                        s_A = model.sigma.a0;
                    % Pseudo-observation noise (anisotropic, inhomogeneous in space)
                    case 'SVDfull'
                        % Compute the basis and variance
                        [tmp_U, tmp_S, tmp_A] = model.sigma.noise_generator.vector_noise_basis_2d(...
                            s_w, model.sigma.N_obs, model.sigma.scaling)
                        % At that point A is the *variance* of small scale
                        % velocity, from which we compute characteristic time.
                        tmp_tau = characteristic_time(model, tmp_A);
                        % Scale A:
                        % we multiply it by tau, and we're done with A.
                        s_A = bsxfun(@times, tmp_A, tmp_tau);
                        % Scale tmp_S:
                        %   i) 1/sqrt(dt) scaling of tmp_S is due to implementation
                        % specificities, see also further down section "Transport of tracer".
                        % It should not be applied to tmp_A, only to tmp_S,
                        % therefore it is not included in svd_noise_basis()
                        % scaling factor. Indeed:
                        % - the BM has a sqrt(dt) factor;
                        % - the advection terms further down have dt in factor,
                        %   so that sigma_dBt neds a (1/dt) factor to
                        %   compensate: we compute sigma_dBt_over_dt;
                        % The combined factor is then sqrt(dt)/dt = 1/sqrt(dt).
                        %   ii) there is also the sqrt(tau) factor. We already
                        % applied a tau factor to A, this is the corresponding one for S.
                        %   So, the overall factor is sqrt(tau/dt), which we 
                        % apply to S (instead of sigma_dBt_over_dt) to minimize computations.
                        tmp_S = tmp_S.*sqrt(tmp_tau/model.advection.dt_adv);
                        % Rebuild the noise sigma.dBt.(1/dt) 
                        s_sigma_dBt_over_dt = model.sigma.noise_generator.draw_noise(...
                            tmp_U, tmp_S);
                    % Noise compensating the hyper-viscoty term
                    case 'Hypervis' %[WIP]   
                        % compute the stochastic advection term from
                        % current buoyancy
                        tmp_Ab = real(ifft2( model.sigma.fft_StoAdv .* s_fft_b )); 
                        % Scramble to generate ensemble of pseudo-observations
                        tmp_Ab_obs = model.sigma.scrambler.scramble_scalar(tmp_Ab, ...
                                                                           model.sigma.N_obs, ...
                                                                           true);
                        % dealias
                        tmp_fft_Ab_obs = fft2( tmp_Ab_obs );
                        tmp_Ab_obs = real(ifft2( bsxfun(@times, tmp_fft_Ab_obs, model.grid.k_aa.mask) ));
                        tmp_Ab_obs = reshape(tmp_Ab_obs, [prod(model.grid.MX), model.sigma.N_obs])
                        % remove mean
                        tmp_Ab_obs = bsxfun(@minus, tmp_Ab_obs, mean(tmp_Ab_obs,2));
                        % apply SVD
                        [tmp_Phi, tmp_lambda, ~] = svd(tmp_Ab_obs, 'econ');
                        tmp_lambda = diag(tmp_lambda)./sqrt(double(model.sigma.N_obs) - 1.); %rescale to get proper variance
                        % draw noise
                        tmp_gamma = bsxfun(@times, randn(numel(tmp_lambda), 1), tmp_lambda); 
                        % so this is the (sigma_dBt. grad(b)) term
                        % (stochastic advection)
                        s_sigma_dBt_gradb = sqrt(model.advection.dt_adv).*reshape(tmp_Phi*tmp_gamma, model.grid.MX);
                        s_fft_sigma_dBt_gradb = fft2( s_sigma_dBt_gradb ); %we work in Fourrier
                        s_fft_sigma_dBt_gradb(model.grid.k.ZM(1),:) = 0.; %so we dealias the single high freq.
                        s_fft_sigma_dBt_gradb(:,model.grid.k.ZM(2)) = 0.;
                        s_sigma_dBt_over_dt = 0.; %we directly compute the advection term above, so the noise velocity is zero
                        s_A = 0.; % because the stochastic diffusion is, by construction, the HV term
                    case 'WavHypervis' %[WIP]
                        % compute the stochastic advection term from
                        % current buoyancy
                        tmp_Ab = real(ifft2( model.sigma.fft_StoAdv .* s_fft_b ));
                        % apply 2D wavelet decomposition
                        tmp_wAb = FWT2_PO(tmp_Ab, model.sigma.Jmin, ...
                                         model.sigma.qmf);
                        % multiply each coefficient by random normal variable
                        tmp_wAb = tmp_wAb.*randn(model.grid.MX);
                        % and rebuild the noise term
                        s_sigma_dBt_gradb = IWT2_PO(tmp_wAb, model.sigma.Jmin, ...
                                                    model.sigma.qmf);
                        s_fft_sigma_dBt_gradb = fft2( s_sigma_dBt_gradb ); %we work in Fourrier
                        s_fft_sigma_dBt_gradb(model.grid.k.ZM(1),:) = 0.; %so we dealias the single high freq.
                        s_fft_sigma_dBt_gradb(:,model.grid.k.ZM(2)) = 0.;
                        s_sigma_dBt_over_dt = 0.; %we directly compute the advection term above, so the noise velocity is zero
                        s_A = 0.; % because the stochastic diffusion is, by construction, the HV term                                          
                    % No noise?! Just in case...
                    case {'None', ''}
                        s_sigma_dBt_over_dt = 0.;
                    % Error                        
                    otherwise
                        error('SQGMU:fft_advection_sto:ValueError', ...
                              '"%s" is not a valid noise model', ...
                              model.sigma.type_noise);
                end
            else % Deterministic case
                s_sigma_dBt_over_dt = 0.;
                s_A = 0.; %[WIP] s_A>0 enable laplacian diffusion on deterministic ensemble
            end
            
            %% Adding time-correlated and time decorrelated velocity
            s_wf = s_w + s_sigma_dBt_over_dt; % wf is the "full" velocity
            
            %% Transport of tracer
            if model.is_stochastic
                % Euler scheme
                % Note: dt in factor of deriv_fft_advection()requires 
                % appropriate compensation on sigma_dBt above.
                s_fft_b = s_fft_b ...
                    + model.advection.dt_adv * ...
                      deriv_fft_advection(model, s_fft_b, s_wf, s_A) ...
                    - s_fft_sigma_dBt_gradb;
            else
                % Runge-Kutta 4 scheme
                s_fft_b = RK4_fft_advection(model, s_fft_b, s_wf, s_A);
            end
            
            %% Last step? save variables for this sample
            if tt==end_t
                % common variables
                fft_buoy_part(:,:,sampl) = s_fft_b; % buoyancy
                w(:,:,:,sampl) = s_w; % large scale (time-correlated) velocity
                if has_probes
                    probes{sampl} = s_probe; % [WIP] probe
                end
                % model-dependent variables
                if model.is_stochastic
                    sigma_dBt_over_dt(:,:,:,sampl) = s_sigma_dBt_over_dt; % noise
                    switch model.sigma.type_noise
                        case {'SVDfull'}
                            A(:,:,:,sampl) = s_A; % diffusion tensor
                        case {'Hypervis', 'WavHypervis'}
                            sigma_dBt_gradb(:,:,sampl) = s_sigma_dBt_gradb;
                    end
                end
            end
        end
    end
    % set current GLOBAL time (end of the interval)
    t = end_t;
    
    %% Discard particles which have blown up
    iii = isnan(fft_buoy_part) | isinf(abs(fft_buoy_part));
    if any(iii(:))
        iii=any(any(any(iii,3),2),1);
        if all(iii(:))
            error('the simulation has blown up');
        end
        nb_dead_pcl = sum(iii);
        warning([ num2str(nb_dead_pcl) ' particle(s) on ' num2str(N_ech) ...
            ' have(s) blown up and' ...
            ' are(is) resampled uniformly on the set of the others particles']);
        N_ech_temp = N_ech - nb_dead_pcl;
        fft_buoy_part(:,:,:,iii)=[];
        iii_sample = randi(N_ech_temp,nb_dead_pcl,1);
        for k=1:nb_dead_pcl
            fft_buoy_part(:,:,:,end+1) = fft_buoy_part(:,:,:,iii_sample(k));
        end
    end
    clear iii
    
    %% Output data
    output()
    
    %% move on
    t = t + 1;
end
if model.verbose
    fprintf(1, '%s - end of simulation\n', datestr(datetime));
end

%% [WIP] save the probe
if has_probes
    output_file = fullfile(model.output.folder_simu, 'files', 'probe');
    save(output_file, 'model', 'probes');
    if model.verbose
            fprintf(1, '\tsaved probes archive: %s\n', output_file);
    end
end

%% nested output function to it can access variables right away.
function output()
    %% function output()
    %
    % Created by P. DERIAN 2016-10-25

    % Time information
    nhours = t*model.advection.dt_adv/3600.; % number of hours (float) 
    ndays = nhours/24.; % Number of days (float)
    hour_str = num2str(floor(nhours)); % integer number here
    day_str = num2str(floor(ndays)); % integer number here
    file_index = day_str; %[TODO] choose the best? how? 
    if model.verbose
        fprintf(1, '%s - iter #%d output: %.3f day(s) of advection\n', ...
                datestr(datetime), t, ndays);
    end
    % Plots
    if model.output.plot_results
        fct_plot(model,fft_buoy_part,file_index);
    end
    % Save file
    output_file = fullfile(model.output.folder_simu, 'files', file_index);
    % first generic variables
    rstate = rng; % retrieve the random generator data
    save(output_file, 'model','t','fft_buoy_part','w', 'rstate');
    % then model-specific variables
    if model.is_stochastic
        switch model.sigma.type_noise
            case 'Spectrum'
                save(output_file, 'sigma_dBt_over_dt', '-append');
            case {'SVD', 'SVDfull', 'Learning_SVD'}
                save(output_file, 'sigma_dBt_over_dt', 'A', '-append');
            case {'Hypervis', 'WavHypervis'}
                save(output_file, 'sigma_dBt_gradb', '-append');
        end
    end
    % Print message
    if model.verbose
        fprintf(1, '\tsaved mat archive: %s\n', output_file);
    end
end

end

