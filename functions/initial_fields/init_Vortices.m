function b_S = init_Vortices(model,X,Y)
% Generate a 2D buoyancy field with two cold cyclones and two warm
% anticyclones
% 

Lx = model.grid.dX(1) * model.grid.MX(1);
ee=4;
odg_b = model.odg_b;
x=X(:,1);
y=Y(1,:);
sigma= 2 * Lx/15; % Length scale close to the Rossby radius

% Warm anticyclones   
center1x=x(1/4*model.grid.MX(1)+1);
center1y=y(1/4*model.grid.MX(2)+1);
b_S = + exp( -1/(2*sigma^2)* (ee*(X-center1x).^2+(Y-center1y).^2) );
center2x=x(3/4*model.grid.MX(1)+1);
center2y=y(1/4*model.grid.MX(2)+1);
b_S = b_S ...
         + exp( -1/(2*sigma^2)* (ee*(X-center2x).^2+(Y-center2y).^2) );

% Cold cyclones     
center1x=x(1/4*model.grid.MX(1)+1);
center1y=y(3/4*model.grid.MX(2)+1);
b_S = b_S ...
         - exp( -1/(2*sigma^2)* (ee*(X-center1x).^2+(Y-center1y).^2) );
center2x=x(3/4*model.grid.MX(1)+1);
center2y=y(3/4*model.grid.MX(2)+1);
b_S =  b_S ...
         - exp( -1/(2*sigma^2)* (ee*(X-center2x).^2+(Y-center2y).^2) );

% Specify the amplitude of the buoyancy
b_S = odg_b* b_S;

