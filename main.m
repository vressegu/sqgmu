%%%%%%%%%%%%%%%%%%%%
%%% Main 
%%%%%%%%%%%%%%%%%%%%
% this initialization is skipped in compile mode
if ~isdeployed
    init;
end

% set/load model
model = set_model();

% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model);

% Advection
[fft_buoy_final, model] = fct_fft_advection_sto(model, fft_buoy);
