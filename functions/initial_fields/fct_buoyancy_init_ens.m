function [fft_sst,model] = fct_buoyancy_init_ens(model)
%% [fft_sst,model] = fct_buoyancy_init_ens(model)
% Create an (ensemble of) initial buoyancy field
%
% Modified by P. DERIAN 2016-10-12: resolution is a member of model.
% Modified by P. DERIAN 2017-04-10: added 'Scramble', 'ScrambleSc' randomizations.
% Modified by P. DERIAN 2017-05-02: 
%   added 'NoiseSVD' enesemble randomization; 
%   added 'type_rand' parameter to discriminate between initial condition
%       and randomization.;
%   added 'Spectrum' ensemble randomization;

%% Grid
n=model.resolution;
m=n;
Lx=1e6;
dx=Lx/n;
dy=dx;
x= dx*(0:n-1);
y= dy*(0:m-1);
model.grid.origin=[0 0];
model.grid.x_ref=x;
model.grid.y_ref=y;
[x,y]=ndgrid(x,y);
model.grid.dX=[dx dy];
MX=[n m];
model.grid.MX=MX;

%% Spatial buoyancy field
% Switch according to the randomization process
switch model.init.type_rand
    
    % No randomization
    %-----------------
    case {'None', ''}
        % Switch according to the initial configuration
        switch model.init.type_data
            case 'Vortices'
                b_S = init_Vortices(model,x,y);
                model.init.label_data = 'Vortices';
            case 'Vortices2' %well-periodized vortices
                b_S = init_Vortices2(model,x,y);
                model.init.label_data = 'Vortices2';
            case 'Perturbed_vortices'
                b_S = init_Perturbed_vortices(model,x,y);
                model.init.label_data = 'Perturbed_vortices';
            case 'PerturbVortices2' %[WIP]
                 b_S = init_PerturbVortices2(model,x,y);
                 model.init.label_data = 'PerturbVortices2';
            case 'Spectrum'
                b_S = init_Spectrum(model);
                model.init.label_data = 'Spectrum';
            otherwise
                error('sqgmu:fct_buoyancy_init:InvalidParameter', ...
                      'The type of initial condition "%s" is unknown', model.init.type_data);
        end
        % replicate
        b_S = repmat(b_S, [1,1,1,model.advection.N_ech]);
        % Issue warning if N_ech>1 and not stochastic, as all particles
        % will be strictly identical
        if (~model.is_stochastic) && (model.advection.N_ech>1)
            warning('sqgmu:fct_buoyancy_init:DeterministicEnsemble',...
                    'The simulation is deterministic (is_stochastic=%d), there is more than 1 realization (N_ech=%d) but no initial randomization (type_rand=%s): all realizations will be strictly identical', ...
                    model.is_stochastic, model.advection.N_ech, model.init.type_rand);
        end
        
    % Scramble from high-res observation of the initial condition
    %-----------------
    case 'Scramble'
        b_S = ens_Scramble(model);
        model.init.label_data = sprintf('%s_ScrambleInit', model.init.type_data);
    
    % Same as Scramble but scaling the variance manually
    %-----------------
    case 'ScrambleSc'
        b_S = ens_ScrambleSc(model);
        model.init.label_data = sprintf('%s_ScrambleScInit', model.init.type_data);         
    
    % Estimate the noise structure from SVD
    %-----------------
    case 'SVDnoise'
        b_S = ens_SVDnoise(model);
        model.init.label_data = sprintf('%s_SVDnoiseInit', model.init.type_data); 
    
    % Add a random spectrum
    %-----------------
    case 'Spectrum'
        b_S = ens_Spectrum(model);
        model.init.label_data = sprintf('%s_SpectrumInit', model.init.type_data);         
   
    % Undefined
    %-----------------
    otherwise
        error('sqgmu:fct_buoyancy_init:InvalidParameter', ...
              'The type of initial randomization "%s" is unknown', model.init.type_rand);
end

%% Fourier transform of the buoyancy field
fft_sst = fft2(b_S);

