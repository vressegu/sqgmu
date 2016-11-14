function b_S = init_Vortices2(model,X,Y)
%% Generate a 2D buoyancy field with two cold cyclones and two warm
% anticyclones
% 
% Note: this is similar to init_Vortices(), but this one is properly
% periodized to avoid gradient errors near domain boundaries.
%
% Written by P. DERIAN 2016-10-12.

% domain
Lx = model.grid.dX(1) * model.grid.MX(1);
ee=4;
odg_b = model.odg_b;
x=X(:,1);
y=Y(1,:);
sigma= 2 * Lx/15; % Length scale close to the Rossby radius

b_S = zeros(model.grid.MX);

%figure; hold on;
for px=[-1 0 1]
    for py = [-1 0 1]
        tmp_X = X + px*Lx;
        tmp_Y = Y + py*Lx;
        % Warm anticyclones   
        center1x=x(1/4*model.grid.MX(1)+1);
        center1y=y(1/4*model.grid.MX(2)+1);
        tmp_b_S = + exp( -1/(2*sigma^2)* (ee*(tmp_X-center1x).^2+(tmp_Y-center1y).^2) );
        center2x=x(3/4*model.grid.MX(1)+1);
        center2y=y(1/4*model.grid.MX(2)+1);
        tmp_b_S = tmp_b_S ...
                 + exp( -1/(2*sigma^2)* (ee*(tmp_X-center2x).^2+(tmp_Y-center2y).^2) );
        % Cold cyclones     
        center1x=x(1/4*model.grid.MX(1)+1);
        center1y=y(3/4*model.grid.MX(2)+1);
        tmp_b_S = tmp_b_S ...
                 - exp( -1/(2*sigma^2)* (ee*(tmp_X-center1x).^2+(tmp_Y-center1y).^2) );
        center2x=x(3/4*model.grid.MX(1)+1);
        center2y=y(3/4*model.grid.MX(2)+1);
        tmp_b_S =  tmp_b_S ...
                 - exp( -1/(2*sigma^2)* (ee*(tmp_X-center2x).^2+(tmp_Y-center2y).^2) );                
        b_S = b_S + tmp_b_S;
        %imagesc(tmp_X(:,1), tmp_Y(1,:), tmp_b_S');
    end
end

% Specify the amplitude of the buoyancy
b_S = odg_b* b_S;

