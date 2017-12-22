function b_S = init_Wei(model,X,Y)
% Generate a 2D buoyancy field similiar to the simulation of Wei Pan
% 

Lx = model.grid.dX(1) * model.grid.MX(1);
odg_b = model.odg_b;
x=X(:,1)/Lx;
y=Y(1,:)/Lx;

b_S = sin(8.0*pi*x)*sin(8.0*pi*y) ...
    + 0.4*cos(6.0*pi*x)*cos(6.0*pi*y) ...
    + 0.3*cos(10.0*pi*x)*cos(4.0*pi*y)...
    + 0.02*sin(2.0*pi*x) + 0.02*sin(2.0*pi*y);

% Specify the amplitude of the buoyancy
b_S = odg_b* b_S;

