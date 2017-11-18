function b_S = init_Constantin_case1(model, X, Y, varargin)
%% Generate a 2D buoyancy field with two cold cyclones and two warm anticyclones
% From Conantiin and al 1994
%

% domain
Lx = model.grid.dX(1) * model.grid.MX(1);
X=2*pi*X/Lx;
Y=2*pi*Y/Lx;

b_S = sin(X).*sin(Y) + cos(Y);

% Specify the amplitude of the buoyancy
odg_b = model.odg_b;
b_S = odg_b* b_S;

end
