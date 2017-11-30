function dt = fct_CFL(model,w)
% Compute the CFL
%


% CFL of the diffusion (or CFL of the white noise advection)
dX2=(model.grid.dX /pi).^2;
bound1=2/model.sigma.a0*prod(dX2)/sum(dX2);

% CFL of the (large-scale) advection
dX=permute(model.grid.dX,[1 3 2]);
bound2=sum(bsxfun(@times,abs(w),pi./dX),3);
bound2=max(bound2(:));
bound2=1/bound2/4;

% CFL of the hyperviscosity
bound3=1/model.advection.HV.maxVal*(prod(dX2)/sum(dX2)) ^ ...
    (model.advection.HV.order/2);
clear dX dX2

% Minimum of the CFL
dt = min([bound1 bound2 bound3]);
clear bound1 bound2 bound3
if model.sigma.sto 
    dt=dt/2;
    % Further constraint on dt due to the use of a (simple) Euler scheme
    % for the SPDE
end
if model.advection.Smag.spatial_scheme
    dt=dt/2;
end

% if model.plots
%     warning('BIDOUILLE SUR DT');
%     dt = dt/10;
%     %     dt = dt/ (20 * model.sigma.Smag.kappamax_on_kappad^2 );
% end

model.advection.dt_adv = dt;

%dt