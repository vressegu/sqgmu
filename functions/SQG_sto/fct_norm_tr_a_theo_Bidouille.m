function norm_tr_a_theo = fct_norm_tr_a_theo_Bidouille(model,kappa1,kappa2)
% Compute the theoretical trace of the variance tensor a
%

norm_tr_a_theo = 1/10 *...
    fct_norm_tr_a_theo_Low_Pass_w_Slope(model,kappa1,kappa2,3/2);

end
