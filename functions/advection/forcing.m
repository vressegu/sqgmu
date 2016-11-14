function F = forcing(model)
%FORCING computes the forcing term for the large-scale velocity
%   Detailed explanation goes here
%
% Written by P. DERIAN 2016-10-06.

switch model.advection.forcing.type_forcing
    case 'Jet'
        % parameters
        y = model.grid.y;
        y0 = y(model.grid.MX(2)/2); % jet center
        std = model.advection.forcing.std;
        amplitude = model.advection.forcing.amplitude;
        % compute the profile
        Fx = amplitude.*exp(-(y-y0).^2./std^2);
        % replicate for all x
        F = zeros(model.grid.MX(1), model.grid.MX(2), 2);
        F(:,:,1) = repmat(Fx, model.grid.MX(1), 1);
    case 'None'
        F = 0.;
    otherwise
        error('SQGMU:forcing:ValueError', '"%s" is not a valid forcing', model.advection.forcing.type_forcing);
end
end

