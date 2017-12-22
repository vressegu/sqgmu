function [fft_sst,model] = fct_buoyancy_init(model,resolution)
% Create an initial buoyancy field
%

%% Grid
n=resolution;
m=n;
Lx=1e6;
% switch dynamics
%     case 'SQG'
%         Lx=1e6;
%     case '2D'
%         Lx=1;        
%     otherwise
%         error('Unknown type of dynamics');
% end
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
switch model.type_data
    case 'Zero'
        b_S = zeros(MX);
    case 'Vortices'
        b_S = init_Vortices(model,x,y);
    case 'disym_Vortices'
        b_S = init_disym_Vortices(model,x,y);
    case 'Vortices2' %well-periodized vortices
        b_S = init_Vortices2(model,x,y);
    case 'Perturbed_vortices'
        b_S = init_Perturbed_vortices(model,x,y);
    case 'PerturbVortices2' %[WIP]
        b_S = init_PerturbVortices2(model,x,y);
    case 'Constantin_case1' 
        b_S = init_Constantin_case1(model,x,y);
    case 'Constantin_case2' 
        b_S = init_Constantin_case2(model,x,y);
    case 'Spectrum'
        b_S = init_Spectrum(model);
    otherwise
        error('this type of initial condition is unknown')
end

%% Fourier transform of the buoyancy field
fft_sst = fft2(b_S);

