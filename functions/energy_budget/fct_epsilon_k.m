function epsi2 = fct_epsilon_k(model,fft_b,day)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end


%% Grid of wave vectors
% kx = model.grid.k.kx; %"normal" grid
% ky = model.grid.k.ky;
%k2 = model.grid.k.k2;
ZM = model.grid.k.ZM; %index of modes ot be zero'ed out
ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
iky_aa = model.grid.k_aa.iky;
%k2_aa = model.grid.k_aa.k2;
mask_aa = model.grid.k_aa.mask; %anti-aliasing mask
mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask
% % if model.advection.Smag.bool
% % %     ikx_aa_LS = model.grid.k_aa_LS.ikx; %anti-aliased grid for gradient
% % %     iky_aa_LS = model.grid.k_aa_LS.iky;
% % %     k2_aa_LS = model.grid.k_aa_LS.k2;
% %     mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask
% % end
% ikx = 1i*model.grid.k.kx;
% iky = 1i*model.grid.k.ky;


M_kappa=min(model.grid.MX);
P_kappa= M_kappa/2;
d_kappa = 2*pi/sqrt(prod(model.grid.MX.* model.grid.dX));
kappa= d_kappa * ( 0:(P_kappa-1) ) ;
% M_kappa=min(model.grid.MX);
% P_kappa= M_kappa/2;
% kappa=1/(M_kappa)* (0:(P_kappa-1)) ;
% kappa=2*pi*max(1./model.grid.dX)*kappa;
k = model.grid.k.k;
kx_plot = k(:,1); kx_plot(ZM(1)) = d_kappa * P_kappa ;

%% Filtering of intial buoyancy
if model.sigma.hetero_energy_flux_prefilter
    fft_b = bsxfun(@times, ...
        model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
        fft_b);
end

%% Time-correlated velocity
fft_w = SQG_large_UQ(model, fft_b);
w=real(ifft2(fft_w));

% dealisasing of the velocity: FT, apply mask, iFT
ft_w = fft2( w );
w_aa = real(ifft2( bsxfun(@times, ft_w, mask_aa) ));
    

%% Deterministic case or other model than SelfSim_from_LS
if ~ model.sigma.sto | ...
        ~ strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
    %model.sigma.km_LS = kappa(2);
    [~,~,~,....
        ~,~,...
        model.sigma.km_LS] = ...
        fct_sigma_spectrum_abs_diff(model,fft_w,false);
end
if ~ model.sigma.sto
    model.sigma.kappamin_on_kappamax = 1/2;
    model.sigma.kappaLS_on_kappamax = 1/8;
end
%% Gradient of b, anti-aliased
% in Fourier space, de-aliased, then in physical space.
fft_gradb_aa(:,:,1,:) = ikx_aa.*fft_b;
fft_gradb_aa(:,:,2,:) = iky_aa.*fft_b;
gradb_aa = real(ifft2(fft_gradb_aa));

%%
epsilon = nan([P_kappa 1 4]);
epsilon_smoothK = nan([P_kappa 1 4]);
PI_loc_sum = zeros([model.grid.MX 4]);
nb_term_PI_loc_sum = 0;
for i_kappa = 1:P_kappa
    % Local wave number
    kappa_local = kappa(i_kappa);
    
    %% Separating large scales and small-scales
%     steep_filter = (k <= kappa_local);
%     steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
%     steep_filter(:,ZM(2)) = 0.;
%     figure(41);
%     subplot(1,2,1);plot(kx_plot,steep_filter(:,1),'b');
%     alpha = 1.;
%     %alpha = 36.;
%     order = 19.;
% %     %ratio_mask_LS = 1;
% %     steep_filter2 = exp(-alpha*( 2./kappa_local ... 
%     steep_filter2 = exp(-alpha*( 1/(eps+kappa_local) ... 
%      	 .* k ).^order );
%     steep_filter2(ZM(1),:) = 0.; %de-alias the single high freq
%     steep_filter2(:,ZM(2)) = 0.;
% %     steep_filter2 = exp(-alpha*( (2./(2*pi/sqrt(prod(model.grid.dX)))) ... 
% %      / ratio_mask_LS .* k ).^order );
%     subplot(1,2,2);plot(kx_plot,steep_filter2(:,1),'r');
%     drawnow;
%     pause(0.2)
    
    alpha = 1.;
    order = 19.;
    steep_filter = exp(-alpha*( 1/(eps+kappa_local) ... 
     	 .* k ).^order );
    steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
    steep_filter(:,ZM(2)) = 0.;
    
    b_LS = real(ifft2( steep_filter .* fft_b ));
%     fft_gradb_aa_SS = bsxfun( @times, (1 - steep_filter) , fft_gradb_aa ) ;
%     fft_gradb_aa_SS(ZM(1),:,:,:) = 0.; %de-alias the single high freq
%     fft_gradb_aa_SS(:,ZM(2),:,:) = 0.;
%     gradb_aa_SS =  real(ifft2(fft_gradb_aa_SS ));
    gradb_aa_SS =  real(ifft2(fft_gradb_aa ));
    
    %% Advection term
    
    % Advective term in physical space
    wgradT(:,:,1) = sum(bsxfun(@times,w_aa,gradb_aa_SS),3);
    wgradT(:,:,2) = real(ifft2( steep_filter .* fft2(wgradT(:,:,1)) ));
    wgradT(:,:,3) = sum(bsxfun(@times,w_aa,gradb_aa),3);
    wgradT(:,:,4) = real(ifft2( steep_filter .* fft2(wgradT(:,:,3)) ));
    
    % Spatially local flux at scale kappa
    PI_loc = bsxfun( @times, b_LS , wgradT);
    
    PI_loc = fft2(PI_loc);
    PI_loc(ZM(1),:,:,:) = 0.; % de-alias the single high freq
    PI_loc(:,ZM(2),:,:) = 0.; % de-alias the single high freq
    PI_loc = real(ifft2(PI_loc));
    
    % Averaged Spatially local flux at scale kappa
    if  kappa_local > model.sigma.km_LS ...
             &  kappa_local <= ...
             (model.sigma.kappaLSforEspi_on_kappamin * model.sigma.kappamin_on_kappamax * kappa(end) - d_kappa)
%              ( model.sigma.kappamin_on_kappamax * kappa(end) - d_kappa)
% %              &  kappa_local < model.sigma.kappamin_on_kappamax * kappa(end)
% %      %        &  kappa_local < model.sigma.kappaLS_on_kappamax * kappa(end)
        nb_term_PI_loc_sum = nb_term_PI_loc_sum + 1;
        PI_loc_sum = PI_loc_sum + PI_loc;
        % PI_loc_smoothK = PI_loc_sum / nb_term_PI_loc_sum;
    %else
    elseif  kappa_local <= model.sigma.km_LS
        PI_loc_smoothK = nan([1 1 4]);
    end
    
    
%     close(figure(18));figure(18);
%     subplot(1,2,1);
%     imagesc(model.grid.x,model.grid.y,PI_loc');
%     axis xy;axis equal; colorbar;
%     PI_loc_LS = real(ifft2(fft2(PI_loc).*model.grid.k_aa.mask));
%     subplot(1,2,2);
%     imagesc(model.grid.x,model.grid.y,PI_loc_LS');
%     axis xy;axis equal; colorbar;
%     drawnow;
%     %
%     %close(figure(19));
%     figure(19);
%     subplot(2,2,1);
%     imagesc(model.grid.x,model.grid.y,real(ifft2(fft_b))');
%     axis xy;axis equal; colorbar;
%     subplot(2,2,2);
%     imagesc(model.grid.x,model.grid.y,b_LS');
%     axis xy;axis equal; colorbar;
%     subplot(2,2,3);
%     imagesc(model.grid.x,model.grid.y,sqrt(sum(gradb_aa_SS.^2,3))');
%     axis xy;axis equal; colorbar;
%     subplot(2,2,4);
%     imagesc(model.grid.x,model.grid.y,wgradT');
%     axis xy;axis equal; colorbar;
%     drawnow;
%     mean(mean(PI_loc,2),1)
%     pause(0.2)
    
    % Global flux at scale kappa
    epsilon(i_kappa,:,:) = mean(mean(PI_loc,2),1);
    epsilon_smoothK(i_kappa,:,:) = mean(mean(PI_loc_smoothK,2),1);

end

PI_loc_smoothK = PI_loc_sum / nb_term_PI_loc_sum;
epsilon_smoothK = bsxfun(@times,mean(mean(PI_loc_smoothK,2),1),ones(P_kappa,1));
%%



% Local wave number
kappa_local =  1/4 * ( model.sigma.km_LS + d_kappa ) + ...
    3/4 * model.sigma.kappamin_on_kappamax * kappa(end) ;
% kappa_local =  (model.sigma.km_LS + ...
%     model.sigma.kappamin_on_kappamax * kappa(end) ) /2;

%% Separating large scales and small-scales

alpha = 1.;
order = 19.;
steep_filter = exp(-alpha*( 1/(eps+kappa_local) ...
    .* k ).^order );
steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
steep_filter(:,ZM(2)) = 0.;

b_LS = real(ifft2( steep_filter .* fft_b ));
fft_gradb_aa_SS = bsxfun( @times, (1 - steep_filter) , fft_gradb_aa ) ;
fft_gradb_aa_SS(ZM(1),:,:,:) = 0.; %de-alias the single high freq
fft_gradb_aa_SS(:,ZM(2),:,:) = 0.;
gradb_aa_SS =  real(ifft2(fft_gradb_aa_SS ));

%% Advection term

% Advective term in physical space
wgradTtemp = sum(bsxfun(@times,w_aa,gradb_aa_SS),3);
wgradT(:,:,5) = real(ifft2( steep_filter .* fft2(wgradTtemp) ));

% Spatially local flux at scale kappa
PI_loc(:,:,5) = bsxfun( @times, b_LS , wgradT(:,:,5));
PI_loc_smoothK(:,:,5) = PI_loc(:,:,5) ;

% % Averaged Spatially local flux at scale kappa
% if  kappa_local >= model.sigma.km_LS ...
%         &  kappa_local < model.sigma.kappamin_on_kappamax * kappa(end)
%     %        &  kappa_local < model.sigma.kappaLS_on_kappamax * kappa(end)
%     nb_term_PI_loc_sum = nb_term_PI_loc_sum + 1;
%     PI_loc_sum = PI_loc_sum + PI_loc;
%     PI_loc_smoothK = PI_loc_sum / nb_term_PI_loc_sum;
%     %else
% elseif  kappa_local < model.sigma.km_LS
%     PI_loc_smoothK = nan([1 1 4]);
% end


% Global flux at scale kappa
epsilon(:,:,5) = mean(mean(PI_loc_smoothK(:,:,5),2),1)*ones(P_kappa,1);
epsilon_smoothK(:,:,5) = epsilon(:,:,5);
    
    %%

% % kappa_VLS = 1/16*kappa(end);
% kappa_VLS = 1/32*kappa(end);

% kappa_VLS = model.sigma.kappa_VLS_on_kappa_LS ...
%     * model.advection.Smag.dealias_ratio_mask_LS * kappa(end);
% % kappa_VLS = 1/4*model.advection.Smag.dealias_ratio_mask_LS*kappa(end);

if model.sigma.hetero_energy_flux_v2
    kappa_VLS = model.sigma.kappa_VLS_on_kappa_LS ...
        * model.sigma.kappamin_on_kappamax * kappa(end);
else
    kappa_VLS = model.sigma.kappa_VLS_on_kappa_LS ...
        * model.advection.Smag.dealias_ratio_mask_LS * kappa(end);
end
alpha = 1.;
order = 19.;
VLS_steep_filter = exp(-alpha*( 1/kappa_VLS ...
    .* k ).^order );
VLS_steep_filter(ZM(1),:) = 0.; %de-alias the single high freq
VLS_steep_filter(:,ZM(2)) = 0.;

%%

if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
        &&    model.sigma.hetero_energy_flux_postfilter
    mask_post_filter = model.grid.k_aa_sigma_hetero_energy_flux_post.mask;
    % mask_post_filter = mask_aa_LS;
    PI_loc_smoothK = real(ifft2( bsxfun(@times, fft2(PI_loc_smoothK), ...
        mask_post_filter) ));
end

%%

%PI_loc_smoothKX = PI_loc_smoothK;

for q=1:5
    PI_loc_smoothK_neg(:,:,q) = - fct_smooth_pos_part(...
        - PI_loc_smoothK(:,:,q));
    PI_loc_smoothK_neg(:,:,q) = ....
        real(ifft2( bsxfun(@times, VLS_steep_filter, ...
        fft2(PI_loc_smoothK_neg(:,:,q)) ) ));
    PI_loc_smoothK_pos(:,:,q) = fct_smooth_pos_part(...
        PI_loc_smoothK(:,:,q));
    PI_loc_smoothK_pre_pos(:,:,q) = PI_loc_smoothK_pos(:,:,q) ...
        + PI_loc_smoothK_neg(:,:,q);
%     PI_loc_smoothKX_pos(:,:,q) = fct_smooth_pos_part(...
%         PI_loc_smoothKX_pre_pos(:,:,q));
end

% Filtering the diffusivity coefficient at large scales
%PI_loc_smoothKX = PI_loc_smoothKX;
PI_loc_smoothKX = fft2(PI_loc_smoothK_pre_pos);
PI_loc_smoothKX = real(ifft2( bsxfun(@times, PI_loc_smoothKX, ...
    mask_aa_LS) ));


for q=1:5
    PI_loc_smoothKX_pos(:,:,q) = fct_smooth_pos_part(...
        PI_loc_smoothKX(:,:,q));
end


%% Coefficient for the Absolute diffusivity by scale
if model.sigma.hetero_energy_flux_averaging_after
    %     PI_loc_smoothKX_pos = PI_loc_smoothKX_pos .^ (1/3) ;
    coef_AbsDif = bsxfun(@times,...
        1./mean(mean(PI_loc_smoothKX_pos.^(1/3),2),1), ...
        PI_loc_smoothKX_pos.^(1/3) );
else
    coef_AbsDif = (bsxfun(@times,...
        1./mean(mean(PI_loc_smoothKX_pos,2),1), ...
        PI_loc_smoothKX_pos )) .^(1/3) ;
end

% coef_AbsDif = bsxfun(@times,...
%     1./mean(mean(PI_loc_smoothKX_pos.^(1/3),2),1), ...
%     PI_loc_smoothKX_pos.^(1/3) );

%% Debug of fct_epsilon_k_onLine

% % %PI_loc_smoothKX_pos_n_2 = ...
% % [~, PI_loc_smoothKX_pos_n_2,PI_loc_smoothKX_2,...
% %     PI_loc_smoothK_neg_2,PI_loc_smoothK_2]...
% %     = ...
%     coef_AbsDiff = ...
%     fct_epsilon_k_onLine(model,fft_b,fft2(w),fft_gradb_aa);
% %     fct_epsilon_k_onLine(model,fft_b,w_aa,fft_gradb_aa);
% 
% PI_loc_smoothKX_2_ref = PI_loc_smoothKX(:,:,2);
% PI_loc_smoothK_neg_2_ref = PI_loc_smoothK_neg(:,:,2);
% PI_loc_smoothK_2_ref = PI_loc_smoothK(:,:,2);
% PI_loc_smoothKX_pos_n_2_ref = PI_loc_smoothKX_pos_n(:,:,2);
% errKX_2 = abs(PI_loc_smoothKX_2 - PI_loc_smoothKX_2_ref) / ...
%     mean(abs(PI_loc_smoothKX_2_ref(:)));
% mean_errorKX_2 = mean(errKX_2(:))
% errK_neg_2 = abs(PI_loc_smoothK_neg_2 - PI_loc_smoothK_neg_2_ref) / ...
%     mean(abs(PI_loc_smoothK_neg_2_ref(:)));
% mean_errorK_neg_2 = mean(errK_neg_2(:))
% errK_2 = abs(PI_loc_smoothK_2 - PI_loc_smoothK_2_ref) / ...
%     mean(abs(PI_loc_smoothK_2_ref(:)));
% mean_errorK_2 = mean(errK_2(:))
% errKX_2_pos_n = abs(PI_loc_smoothKX_pos_n_2 - PI_loc_smoothKX_pos_n_2_ref) / ...
%     mean(abs(PI_loc_smoothKX_pos_n_2_ref(:)));
% mean_errorKX_2_pos_n = mean(errKX_2_pos_n(:))
% 
% errKX_2_pos_n = abs(PI_loc_smoothKX_pos_n_2 - PI_loc_smoothKX_pos_n_2_ref) / ...
%     mean(abs(PI_loc_smoothKX_pos_n_2_ref(:)));
% mean_errorKX_2_pos_n = mean(errKX_2_pos_n(:))

%%


close(figure(18));figure18=figure(18);
width = 20;
height = 18;
set(figure18,'Units','inches', ...
    'Position',[0 0 width height], ...
    'PaperPositionMode','auto');
for q=1:5
    subplot(5,6,1+6*(q-1));
    imagesc(model.grid.x,model.grid.y,PI_loc(:,:,q)');
    axis xy;axis equal; colorbar;
    title('one k');
    subplot(5,6,2+6*(q-1));
    imagesc(model.grid.x,model.grid.y,PI_loc_smoothK(:,:,q)');
    axis xy;axis equal; colorbar;
    title('Smoothink along k');
    subplot(5,6,3+6*(q-1));
    %surf(model.grid.x,model.grid.y,PI_loc_smoothKX_pos');
    imagesc(model.grid.x,model.grid.y,PI_loc_smoothK_pre_pos(:,:,q)');
    axis xy;axis equal; colorbar;
    title('pre >0 part');
    subplot(5,6,4+6*(q-1));
    %surf(model.grid.x,model.grid.y,PI_loc_smoothKX');
    imagesc(model.grid.x,model.grid.y,PI_loc_smoothKX(:,:,q)');
    axis xy;axis equal; colorbar;
    title('Smoothink along (k,x)');
    subplot(5,6,5+6*(q-1));
    %surf(model.grid.x,model.grid.y,PI_loc_smoothKX_pos');
    imagesc(model.grid.x,model.grid.y,PI_loc_smoothKX_pos(:,:,q)');
    axis xy;axis equal; colorbar;
    title('>0 part');
    subplot(5,6,6+6*(q-1));
    imagesc(model.grid.x,model.grid.y,...
        (coef_AbsDif(:,:,q)'));
%        (PI_loc_smoothKX_pos_n(:,:,q)').^(1/3));
    axis xy;axis equal; colorbar;
    title('Coef AbsDif');
end
drawnow;
eval( ['print -depsc ' model.folder.folder_simu ...
    '/Epsilon_k_meth/' day '.eps']);
figure8=figure(8);

% close(figure(18));figure(18);
% subplot(2,2,1);
% imagesc(model.grid.x,model.grid.y,PI_loc_smoothK');
% axis xy;axis equal; colorbar;
% subplot(2,2,2);
% %surf(model.grid.x,model.grid.y,PI_loc_smoothKX');
% imagesc(model.grid.x,model.grid.y,PI_loc_smoothKX');
% axis xy;axis equal; colorbar;
% subplot(2,2,3);
% %surf(model.grid.x,model.grid.y,PI_loc_smoothKX_pos');
% imagesc(model.grid.x,model.grid.y,PI_loc_smoothKX_pos');
% axis xy;axis equal; colorbar;
% subplot(2,2,4);
% imagesc(model.grid.x,model.grid.y,PI_loc_smoothKX_pos_n');
% axis xy;axis equal; colorbar;
% drawnow;

%% Plot
% slope_ref = 0;
% idx_not_inf=~(isinf(log10(epsilon(2:end))) ...
%     | epsilon(2:end)<1e-4*max(epsilon(2:end)) | isinf(kappa(2:end)'));
% idx_not_inf = [false; idx_not_inf];
% % line1= slope_ref * log10(kappa(2:end))  ;
% line1= slope_ref * log10(kappa(2:end))  ;
% offset = -1 + mean(  log10(epsilon(idx_not_inf)')  ...
%     - line1(idx_not_inf(2:end)));
% line1 = line1 + offset;
%ref=10.^line1;
ref=zeros(size(kappa(2:end)));
semilogx(kappa(2:end),ref,'--k');
hold on;
% name_plot =semilogx(kappa(2:end) , epsilon(2:end) ,color);
name_plot =semilogx(kappa(2:end) , epsilon(2:end,:,1) ,'b+-');
name_plot =semilogx(kappa(2:end) , epsilon(2:end,:,2),'b*-');
name_plot =semilogx(kappa(2:end) , epsilon(2:end,:,3) ,'bx-');
name_plot =semilogx(kappa(2:end) , epsilon(2:end,:,4) ,'bs-');
name_plot =semilogx(kappa(2:end) , epsilon(2:end,:,5) ,'bd-');
%semilogx(kappa(2:end) , epsilon_smoothK(2:end) ,'r');
semilogx(kappa(2:end),epsilon_smoothK(end,:,1)*ones(size(kappa(2:end))),'-+r')
epsi1 = epsilon_smoothK(end,:,1)
semilogx(kappa(2:end),epsilon_smoothK(end,:,2)*ones(size(kappa(2:end))),'r-*')
epsi2 = epsilon_smoothK(end,:,2)
semilogx(kappa(2:end),epsilon_smoothK(end,:,3)*ones(size(kappa(2:end))),'r-x')
epsi3 = epsilon_smoothK(end,:,3)
semilogx(kappa(2:end),epsilon_smoothK(end,:,4)*ones(size(kappa(2:end))),'r-s')
epsi4 = epsilon_smoothK(end,:,4)
semilogx(kappa(2:end),epsilon_smoothK(end,:,5)*ones(size(kappa(2:end))),'r-d')
epsi5 = epsilon_smoothK(end,:,5)
ax=axis;
ax(4)=max([max(epsilon(2:end,:,:),[],3); ref']);
% min_ax= 10 ^(slope_ref * log10(kappa(2)*512/2) + offset) ;
% ax(3) = (model.odg_b/(1e-3))^2 * ...
%     6e-2*(kappa(2)/kappa(end))*min([max(epsilon); max(ref)']);
% % ax(3)=6e-2*(kappa(2)/kappa(end))*min([max(epsilon); max(ref)']);
% ax(3) = min( [ax(3) min(ref) min_ax]);
ax(1:2)=kappa(2)*[1 min(model.grid.MX)/2];
if ax(4)>0
    axis(ax)
end
ax=axis;
semilogx((d_kappa+model.sigma.km_LS)*[1 1],...
    [min([ax(3) ref]) ax(4)],'k--')
semilogx(model.sigma.kappamin_on_kappamax * kappa(end)*[1 1],...
    [min([ax(3) ref]) ax(4)],'k-.')
