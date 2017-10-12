function [fft_buoy_part, model] = fct_fft_advection_sto(model,  fft_buoy_part)
% Advection of buoyancy using SQG or SQG_MU model
%

%% Folder to save plots and files
if model.advection.HV.bool
    add_subgrid_deter = '_HV';
elseif model.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
    %add_subgrid_deter = [];
end
if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
        model.advection.Smag.bool
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
end

if isinf(model.sigma.k_c) % Deterministic case
    model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
else % Stochastic case
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_MU' add_subgrid_deter '/' model.type_data ];
end
if model.advection.forcing.bool
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_forced_turb' ];
else
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_free_turb' ];
end
model.folder.folder_simu = [ model.folder.folder_simu ...
    '/' num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];
% Create the folders
fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
N_ech=model.advection.N_ech;

%% Initialisation of the spatial fields

% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_buoy_part);
w=real(ifft2(fft_w));

% figure;imagesc(real(ifft2(fft_buoy_part))');axis xy;axis equal;colorbar;
% figure;imagesc(w(:,:,1)');axis xy;axis equal;colorbar;

% Create several identical realizations of the intial buoyancy
fft_buoy_part = repmat(fft_buoy_part(:,:,1),[1 1 1 model.advection.N_ech]);

%% Choice of the variance tensor a
% Variance tensor
model.sigma.a0 = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
% Diffusion coefficient
model.advection.coef_diff = 1/2 * model.sigma.a0;

%% Forcing
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %     Ly = model.grid.MX(2) * model.grid.dX(2);
    %     [~,Y]=ndgrid(model.grid.x,model.grid.y);
    %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         sin( 2 * pi * model.advection.forcing.freq_f * Y / Ly);
    %     clear Ly X Y
    
    
    model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
    model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
    [X,Y]=ndgrid(model.grid.x,model.grid.y);
    
    
    
    switch model.dynamics
        case 'SQG'
            U_caract = ...
                model.odg_b / model.physical_constant.buoyancy_freq_N;
        case '2D'
            U_caract = ...
                model.odg_b * model.advection.forcing.Ly;
            U_caract = U_caract /20;
            U_caract = U_caract /4;
        otherwise
            error('Unknown type of dynamics');
    end
    %     U_caract = U_caract /10;
    U_caract = U_caract /3;
    
    %     model.advection.U_caract = sqrt(mean(w(:).^2));
    model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
    
    %     model.advection.forcing.on_T = model.advection.forcing.on_T /10;
    
    model.advection.forcing.ampli_forcing = ...
        model.advection.forcing.ampli_forcing * model.odg_b ;
    %     model.advection.forcing.ampli_forcing = ...
    %         model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         model.advection.forcing.on_T;
    
    %     ampli_scale = sqrt( ...
    %         ( sum(model.advection.forcing.freq_f.^2)/2 ) ...
    %             ^(model.sigma.slope_sigma) ...
    %             );
    ampli_scale = 1;
    
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing * ...
                cos( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X ...
                + 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y );
            F = F / model.advection.forcing.on_T;
        case 'Spring'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing ...
                * sin( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X ) ...
                .* sin( 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y );
            %     F = model.advection.forcing.ampli_forcing * ...
            %         sin( 2 * pi * model.advection.forcing.freq_f * Y / ...
            %         model.advection.forcing.Ly);
    end
    
    
    
    model.advection.forcing.F = fft2(F);
    
    if strcmp(model.type_data, 'Zero')
        fft_w = SQG_large_UQ(model,  ...
            model.odg_b / model.advection.forcing.ampli_forcing ...
            * model.advection.forcing.F);
        %         w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
        %         w(:,:,2) = 0;
        %         fft_w = fft2(w);
        if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
            fft_w = fft_w * model.advection.forcing.on_T;
        end
        w = real(ifft2(fft_w));
    end
    
    %
    %     figure;imagesc(model.grid.x,model.grid.y,model.advection.forcing.F');
    %     axis xy; axis equal;
    %
    
    %     model.advection.forcing.F = fft2(F);
    clear Lx Ly X Y on_T U_caract ampli_scale
end

%% Hyperviscosity

% Root mean square Lyapunov
[dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s=(dxUx+dyUy)/2;
d=(dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s d dxUx dxUy dyUx dyUy lambda

% Order of the hyperviscosity
if model.advection.HV.bool
    model.advection.HV.order=8;
else
    model.advection.HV.order=2;
end

% Hyperviscosity coefficient
model.advection.HV.val= ...
    40 * model.advection.lambda_RMS * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;

% if ~model.advection.HV.bool % DNS
if model.advection.Lap_visco.bool &&  strcmp(model.dynamics,'2D')
    if isfield(model.advection, 'forcing') && model.advection.forcing.bool
        model.carct.L_caract = 1/ ...
            sqrt(mean( ( [ model.advection.forcing.Lx ...
            model.advection.forcing.Ly ] ...
            .\ model.advection.forcing.freq_f ).^2));
        warning('Re should be computed with the amplitude of forcing!');
    else
        fft_grad_b = fct_grad(model,fft_buoy_part);
        model.carct.L_caract = sqrt( ...
            sum(abs(fft_buoy_part(:)).^2) / sum(abs(fft_grad_b(:)).^2) );
    end
    model.carct.U_caract = model.odg_b * model.carct.L_caract;
    Re = model.carct.U_caract * model.carct.L_caract ...
        / model.advection.HV.val
    model.carct.Re = Re;
    model.carct.l_Kolmogorov = Re^(-1/2) * model.carct.L_caract;
    if model.carct.l_Kolmogorov < min(model.grid.dX)
        error('The simulation is under resolved');
    end
end

if model.advection.Smag.bool
    model.advection.Smag.coef_Smag = ...
        ( model.advection.Smag.kappamax_on_kappad ...
        * sqrt(prod(model.grid.dX))/pi ) ...
        ^ (5*model.advection.HV.order/4 -1/2);
end

%% Choice of time step : CFL

% CFL of the diffusion (or CFL of the white noise advection)
dX2=(model.grid.dX /pi).^2;
bound1=2/model.sigma.a0*prod(dX2)/sum(dX2);

% CFL of the (large-scale) advection
dX=permute(model.grid.dX,[1 3 2]);
bound2=sum(bsxfun(@times,abs(w),pi./dX),3);
bound2=max(bound2(:));
bound2=1/bound2/4;

% CFL of the hyperviscosity
bound3=1/model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
    (model.advection.HV.order/2);
clear dX dX2

% Minimum of the CFL
dt = min([bound1 bound2 bound3]);
clear bound1 bound2 bound3
if ~ isinf(model.sigma.k_c)
    dt=dt/2;
    % Further constraint on dt due to the use of a (simple) Euler scheme
    % for the SPDE
end

warning('BIDOUILLE SUR DT');
dt = dt/3;

model.advection.dt_adv = dt;

%% Fourier transform of the kernel \tilde sigma
if ~isinf(model.sigma.k_c)
    % Fourier transform of the kernel \tilde sigma up to a multiplicative
    % constant
    [sigma_on_sq_dt, ~, missed_var_small_scale_spectrum ] ...
        = fct_sigma_spectrum(model,fft_w);
    % sigma_on_sq_dt will be used to simulate sigma d B_t
    % missed_var_small_scale_spectrum will be used to set the mulitplicative
    % constant
    
    % Muliplicative constant of the kernel \tilde sigma
    model.sigma.a0_on_dt = model.sigma.a0 / dt;
    sigma_on_sq_dt = sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
        * sigma_on_sq_dt;
    % the factor d=2 is for the dimension d of the space R^d
    % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
    clear missed_var_small_scale_spectrum
    warning('This formula may changed if k inf > la resolution');
else
    sigma_on_sq_dt = 0;
end

%% Possibly reset velocity to zero
% if strcmp(model.type_data, 'Zero')
%     w = zeros(size(w));
% end
clear w fft_w

%% Loop on time

% Used for the first plot
tt_last = -inf;

% Number of time step
N_t = ceil(model.advection.advection_duration/dt);

%% Use a saved files of a former simulation ?
if model.advection.use_save
    day = num2str(model.advection.day_save);
    name_file = [model.folder.folder_simu '/files/' day '.mat']; clear day
    load(name_file)
    
    % Version of matlab
    vers = version;
    year = str2double(vers(end-5:end-2));
    subvers = vers(end-1);
    model.folder.colormap_freeze = ...
        (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
    
    t_ini=t+1;
    tt_last = -inf;
else
    t_ini=1;
end

% Printing some information
fprintf(['The initial condition is ' model.type_data ' \n'])
fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
fprintf(['Time step : ' num2str(dt) ' seconds \n']);
fprintf(['Time of advection : ' num2str(N_t*dt/3600/24) ' days \n']);
fprintf(['Ensemble size : ' num2str(N_ech) ' realizations \n']);

for t=t_ini:N_t
    % for t=1:N_t
    %% Time-uncorrelated velocity (isotropic and homogeneous in space)
    if isinf(model.sigma.k_c) % Deterministic case
        sigma_dBt_dt = 0;
    else % Stochastic case
        % Fourier transform of white noise
        dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
        % Multiplication by the Fourier transform of the kernel \tilde \sigma
        fft_sigma_dBt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
        clear dBt_C_on_sq_dt
        % Homogeneous velocity field
        sigma_dBt_dt = real(ifft2(fft_sigma_dBt));
        clear fft_sigma_dBt
    end
    
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_buoy_part);
    w=real(ifft2(fft_w));
    clear fft_w
    
    %% Adding time-correlated and time decorrelated velocity
    w = w + sigma_dBt_dt;
    %     if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %         w(:,:,1) = w(:,:,1) + Vy;
    %     end
    
    %% Transport of tracer
    if isinf(model.sigma.k_c)
        % Runge-Kutta 4 scheme
        fft_buoy_part = RK4_fft_advection(model,fft_buoy_part, w);
    else
        parfor sampl=1:N_ech
            % Euler scheme
            fft_buoy_part(:,:,:,sampl) = fft_buoy_part(:,:,:,sampl) ...
                + deriv_fft_advection( ...
                model, fft_buoy_part(:,:,:,sampl), w(:,:,:,sampl)) * dt;
        end
    end
    
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
    
    %% Plots and save
    tt = floor(t *dt/ (3600*24)); % Number of days
    if tt > tt_last
        tt_last = tt;
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        day = num2str(floor(t*dt/24/3600));
        
        % Plots
        fct_plot(model,fft_buoy_part,day)
        
        if model.advection.plot_dissip
            fct_plot_dissipation(model,fft_buoy_part,sigma_on_sq_dt,day);
        end
        
        % Save files
        save( [model.folder.folder_simu '/files/' day '.mat'], ...
            'model','t','fft_buoy_part','w','sigma_dBt_dt', ...
            'sigma_on_sq_dt');
    end
end