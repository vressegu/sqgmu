function b_S = init_PerturbVortices2(model, X, Y)

    N = model.init.N_avg;
    % particles weights
    weights = randn(N).^2;
    weights = weights/sum(weights);
    % vortex offsets
    offset = model.init.delta_err*randn(4,2,N);
    offset(:,1,:) = offset(:,1,:)*model.grid.dX(1);
    offset(:,2,:) = offset(:,2,:)*model.grid.dX(2);
    % combine
    b_S = 0.;
    for n=1:N
        b_S = b_S + weights(n)*init_Vortices2(model, X, Y, offset(:,:,n)); 
    end
end