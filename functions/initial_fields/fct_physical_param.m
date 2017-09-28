function model = fct_physical_param(dynamics)
% Specify physical parameters of the algortihm
%

if nargin < 1
    dynamics = 'SQG';
end

% Coriolis parameter
angle_grid = 45/360*2*pi; % rad
OMEGA = 2*pi/(24*60*60);
model.physical_constant.f0 = 2* OMEGA * sin( angle_grid );

% Background density
model.physical_constant.rho=1e3;

% Gravity
model.physical_constant.g=9.81;

switch dynamics
    case 'SQG'        
        % Background stratification
        model.physical_constant.buoyancy_freq_N = ...
            3 * model.physical_constant.f0;
        
        % Amplitude of the buoyancy
        model.odg_b = 1e-3;
        
    case '2D'
        % Amplitude of the vorticity
        model.odg_b = 1e-5;
        
        % Background stratification
        model.physical_constant.buoyancy_freq_N = 1;
        
    otherwise
        error('Unknown type of dynamics');
end
