function [sigma_on_sq_dt,f_sigma,trace_a,....
    slope_w_a_comp_for_estim,mult_offset_spectrum_a_estim,...
    km_LS, ...
    spectrum_a_sigma] = ...
    fct_sigma_spectrum_abs_diff_postprocess(model,ft_w,bool_plot,day)
% - sigma_on_sq_dt is the Fourier transform of the kernel \tilde sigma up to a multiplicative
% constant
% - f_sigma is the Fourier transform of the associted streamfunction
% - trace_a_on_dt measures the total energy of the field and is generally used to
% set the muliplicative constant
% - spectrum_a_sigma is the spectrum
%

%%

LineWidth = 1.3;
MarkerSize = 8;
Color1=[0.8 0.1 0.1];
%             Color1=[0.8 0.0 0.1];
Color2=[0.1 0.0 0.8];
Color3=[0.0 0.5 0.0];
%         Color3=[0.0 0.8 0.2];


%         set(name_plot,'LineWidth',LineWidth,...
%             'MarkerSize',MarkerSize,...
%             'Color',Color2);
%         hold on;

%%

% Average over realizations
ft_w2=mean(abs(ft_w).^2,4);

% Fourier transform norm
ft_w2=sum(ft_w2,3);
% if nargin >2
%     ft2=mean(abs(ft2).^2,4);
%     ft2=sum(ft2,3);
% end

if ~isfield(model.sigma,'slope_sigma_ref')
    warning('No field slope_sigma_ref -> default values');
    switch model.dynamics
        case 'SQG'
            model.sigma.slope_sigma_ref = -5/3;
        case '2D'
            model.sigma.slope_sigma_ref = -3;
        otherwise
            error('Unknown type of dynamics');
    end
end
% slope_ref = model.sigma.slope_sigma;
slope_ref = model.sigma.slope_sigma_ref;
slope_ref_a = (slope_ref-3)/2;

% Get parameters
MX=model.grid.MX;
dX=model.grid.dX;
if any(size(ft_w2)~=MX)
    error('wrong size');
end
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;
% ft_w2(PX(1),:)=0;
% ft_w2(:,PX(2))=0;
ft_w2(PX(1)+1,:,:,:)=0;
ft_w2(:,PX(2)+1,:,:)=0;

%% Wave vector
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
[kx,ky]=ndgrid(kx,ky);
k=sqrt(kx.^2+ky.^2);
k(PX(1)+1,:)=inf;
k(:,PX(2)+1)=inf;
k=k(:);

%% Wave number
M_kappa=min(model.grid.MX);
P_kappa= M_kappa/2;
d_kappa = 2*pi/sqrt(prod(model.grid.MX.* model.grid.dX));
kappa= d_kappa * ( 0:(P_kappa-1) ) ;
% M_kappa=min(MX);
% P_kappa= M_kappa/2;
% %d_kappa = max(1./dX);
% kappa=1/(M_kappa)* (0:(P_kappa-1)) ;
% kappa=2*pi*max(1./dX)*kappa;
% %kappa=2*pi*d_kappa*kappa;

%% Masks associated with the rings of iso wave number
d_kappa = kappa(2) - kappa(1);
idx = sparse( bsxfun(@le,kappa, k ) );
idx = idx & sparse( bsxfun(@lt,k, [ kappa(2:end) kappa(end)+d_kappa ] ) );

kappa = kappa';

%% Spectrum
% Integration over the rings of iso wave number
spectrum_w = idx' * ft_w2(:);
% if nargin>2
%     spectrum_w2 = idx' * ft2(:);
% end

% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% velocity averaged (not just integrated) over the space
spectrum_w = 1/prod(model.grid.MX)^2 * spectrum_w;
% if nargin>2
%     spectrum_w2 = 1/prod(model.grid.MX)^2 * spectrum_w2;
% end

% Division by the wave number step
d_kappa = kappa(2)-kappa(1);
spectrum_w = spectrum_w / d_kappa;
% if nargin>2
%     spectrum_w2 = spectrum_w2 / d_kappa;
% end

% % From buoyancy spectrum to velocity spectrum
% spectrum_w = spectrum_w / model.physical_constant.buoyancy_freq_N^2;

% From velocity spectrum to absolute diffusity by scale
spectrum_w_a = zeros(size(kappa));
spectrum_w_a(2:end) = kappa(2:end).^(-3/2) .* spectrum_w(2:end).^(1/2);

% Time scales
% %tau_k = 1./sqrt( d_kappa * kappa.^2 .* spectrum_w);
% v_k = sqrt( d_kappa * spectrum_w) * sqrt(prod(model.grid.MX));
% tau_k = 1./sqrt( kappa.^3 .* spectrum_w);
% figure(56);loglog(kappa(2:end),v_k(2:end)/3600/24);
% figure(55);loglog(kappa(2:end),tau_k(2:end)/3600/24);

% Compensated absolute diffusity by scale
spectrum_w_a_comp = kappa.^(-slope_ref_a) .* spectrum_w_a;

%% Absosute diffusivity
abs_diff_w = sum(spectrum_w_a) * d_kappa;

% %% Time integral of the dissipation per scale epsilon(k)
% int_epsilon = - cumsum(spectrum) * d_kappa;

%% Estimation of slope of the spectral absolute diffusivity of
% % the large-scale velocity

% Get the large-scale "length scale" km_LS and the spectral window for the
% linear regression
% threshold_k = 1/2;
if ~ model.sigma.sto
    model.sigma.kappamin_on_kappamax = 1/2;
    model.sigma.kappaLS_on_kappamax = 1/8;
    pre_estim_slope=1e-1;
    %     pre_5 = 5e-2;
    model.sigma.kappamin_on_kappamax_estim_slope = ...
        (log(1-pre_estim_slope)/log(pre_estim_slope))...
        ^(2/model.advection.HV.order);
    model.sigma.estim_k_LS = false;
end
threshold_k = model.sigma.kappamin_on_kappamax_estim_slope;
% threshold_k = model.sigma.kappamin_on_kappamax;
threshold_k_LS = model.sigma.kappaLS_on_kappamax;
%threshold_k_LS = 1/8;
spectrum_w_a_comp_cut = spectrum_w_a_comp(kappa < kappa(end)*threshold_k);
% %spectrum_w_a_comp_cut = spectrum_w_a_comp;
% [~,i_first]=max(spectrum_w_a_comp_cut(2:end));
%%
if model.sigma.estim_k_LS
    [~,i_first]=max(spectrum_w_a_comp_cut( [ false; ...
        ( kappa(2:end) < kappa(end)*threshold_k_LS )] ) );
else
    i_first = 2;
    warning('Spectrum slope estimated from the largest scales')
end
%%
%[~,i_first]=max(spectrum_w_a_comp(2:end));
i_first = i_first +1;
iii_k_LS = i_first:length(spectrum_w_a_comp_cut);
spectrum_w_a_comp_for_estim = spectrum_w_a_comp_cut(iii_k_LS);
kappa_cut = kappa(kappa<kappa(end)*threshold_k);
kkk=kappa_cut(iii_k_LS);
km_LS=kkk(1);
offset_w_a = spectrum_w_a(iii_k_LS(1));
offset_w_a_comp = spectrum_w_a_comp(iii_k_LS(1));

% Linear regression
%threshold_k = 1/6;
iii_reliable=~(isinf(abs(spectrum_w_a_comp_for_estim))|...
    isnan(spectrum_w_a_comp_for_estim)|...
    (spectrum_w_a_comp_for_estim/max(spectrum_w_a_comp_for_estim(:)) ...
    <=eps)|...
    (kkk>kappa(end)*threshold_k));
%     (kkk'>kidx(end)/3));
%     (kkk'>10*km));
logspectrum_w_a_comp_for_estim_centered = ...
    log10(spectrum_w_a_comp_for_estim(iii_reliable))...
    - log10(offset_w_a_comp);
logkkk = log10(kkk(iii_reliable)) - log10(km_LS);
% Z = logkkk;
% Z(:,2) = 1;
beta = logkkk \ logspectrum_w_a_comp_for_estim_centered;
beta =min([0 beta]); % Prevent unstable behavior of this parametrisation.
slope_w_a_comp_for_estim = beta + slope_ref_a;
slope_w_comp_for_estim = 2 * slope_w_a_comp_for_estim + 3;
%offset = beta(2);
clear beta


%% 1D Spectrum

% Largest wave number
k_inf = kappa(min(PX));

% Smallest wave number
k0 = model.sigma.kappamin_on_kappamax * k_inf;

% Slope of the absolute diffusitiy distribution by scale
model.sigma.slope_absDif_sigma = (model.sigma.slope_sigma-3)/2;

% Absolute diffusivity by scale with estimated slope
reference_spectrum_a_estim = kappa(2:end) .^ slope_w_a_comp_for_estim ;

% Setting offset
mult_offset_spectrum_a_estim = spectrum_w_a(iii_k_LS(1)) ...
    / reference_spectrum_a_estim(iii_k_LS(1)-1);
%10^offset
% reference_spectrum_a_estim = 10^offset...
%     * reference_spectrum_a_estim;
reference_spectrum_a_estim = mult_offset_spectrum_a_estim...
    * reference_spectrum_a_estim;

% Absolute diffusivity by scale with theoretical slope
reference_spectrum_a = kappa(2:end) .^ slope_ref_a ;
%reference_spectrum_a = kappa(2:end) .^ model.sigma.slope_sigma ;

% Setting offset
mult_offset_spectrum_a = spectrum_w_a(iii_k_LS(1)) ...
    / reference_spectrum_a(iii_k_LS(1)-1) ;
reference_spectrum_a = mult_offset_spectrum_a * reference_spectrum_a;

% % Compute the offset to superimpose the slop -5/3 and the spectrum of w
% idx_not_inf=~(isinf(log10(spectrum_w_a(2:end)))| ...
%     spectrum_w_a(2:end)<1e-4*max(spectrum_w_a(2:end)) | isinf(kappa(2:end)'));
% idx_not_inf= [ false; idx_not_inf ];
% idx_not_inf(end)=false;
% % reference_spectrum = reference_spectrum * 10 .^( ...
% %     mean( log10(spectrum_w(idx_not_inf)') ...
% %     - log10( reference_spectrum(idx_not_inf)) ) );
% mult_offset = 10 .^( ...
%     mean( log10(spectrum_w_a(idx_not_inf)') ...
%     - log10( reference_spectrum_a(idx_not_inf)) ) );
% %reference_spectrum = reference_spectrum * mult_offset;

% Residual absolute diffusivity by scale
%%
if model.sigma.sto
    residual_spectrum_a = ...
        reference_spectrum_a_estim - spectrum_w_a(2:end);
    % residual_spectrum_a = ...
    %     reference_spectrum_a_estim;
    % warning('No residual spectrum')
    %%
    residual_spectrum_a(residual_spectrum_a<0)=0;
    % residual_spectrum_a = ...
    %     reference_spectrum_a_estim ./ spectrum_w_a(2:end);
    % % residual_spectrum_a = ...
    % %     reference_spectrum_a ./ reference_spectrum_a_estim;
    % % % residual_spectrum_a = ...
    % % %     (mult_offset_spectrum_a / mult_offset_spectrum_a_estim) * ...
    % % %     reference_spectrum_a ./ reference_spectrum_a_estim;
    residual_spectrum_a(1:(iii_k_LS(1)-1))=0;
    residual_spectrum_a(1:2)=0;
    
    %% Spectre of random small-scale velocity
    f_sigma = residual_spectrum_a;
    %f_sigma = f_sigma * d_kappa;
    
    % switch model.sigma.type_spectrum
    %     case 'Band_Pass_w_Slope'
    %         f_sigma = reference_spectrum;
    %
    %         % % Parseval ( * prod(model.grid.MX) )
    %         % % and integrated over the space ( * prod(model.grid.MX) )
    %         % f_sigma = prod(model.grid.MX)^2 * f_sigma;
    %
    % Band-pass filter
    idx1 = (kappa(2:end) <= k0);
    idx3 = (kappa(2:end) > k_inf);
    idx2 = ~ (idx1 | idx3);
    
    unit_approx = fct_unity_approx_(sum(idx2));
    % unit_approx = ones([1,sum(idx2)]);
    % warning('Above line modified for debug !!!!!')
    
    f_sigma(idx1 | idx3)=0;
    f_sigma(idx2) = f_sigma(idx2) .* unit_approx';
    
    %     case 'Low_Pass_w_Slope'
    %         f_sigma = (k0^2 + kappa(2:end).^2) .^ (model.sigma.slope_sigma/2) ;
    %     case  'Low_Pass_streamFct_w_Slope'
    %         f_sigma = kappa(2:end).^2 .* ...
    %             (k0^2 + kappa(2:end).^2) .^ (model.sigma.slope_sigma/2-1) ;
    %     case  'BB'
    %         f_sigma = ones(size(kappa(2:end))) ;
    %     case  'Bidouille'
    %         f_sigma = 1/10 * ones(size(kappa(2:end))) ;
    %     otherwise
    %         error('Unknown spectrum type for the small-scale velocity');
    % end
    
    % To remove the 2 pi which appear when we integrate the spectrum over the
    % wave-vector angles and add the (2 pi)^2 which appear when we go from k to
    % 2*pi*k
    f_sigma = f_sigma * (2*pi);
    
    % Division by prod(model.grid.MX) because of the variance white-in-space
    % noise
    f_sigma = 1/prod(model.grid.MX) * f_sigma;
    
    % Multiplication by prod(model.grid.MX) because the spectrum respresent
    % an energy by unit of space and not the sum over the space
    f_sigma = prod(model.grid.MX) * f_sigma;
    
    % Discretisation to go from continuous to discrete Fourier transform
    f_sigma = 1/prod(model.grid.dX) * f_sigma;
    
    % Division by k^2 to get the spectrum of the streamfunction
    f_sigma = f_sigma ./ ( kappa(2:end) .^2 );
    
    % From omnidirectional spectrum to Fourier tranform square modulus
    % Division by k in dimension 2
    f_sigma = f_sigma ./ ( kappa(2:end) );
    
    % % Influence of discretisation
    % f_sigma = f_sigma * d_kappa ;
    % f_sigma = f_sigma / ( prod(MX.*dX) /(2*pi) ) ;
    
    % From square modulus to modulus
    f_sigma = sqrt( f_sigma );
    
    % From 1D function to 2D function
    f_sigma = interp1(kappa,[0; f_sigma],k);
    
    % % Cleaning
    % if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
    %     f_sigma(k<=k0)=0;
    % end
    f_sigma(k>k_inf)=0;
    f_sigma=reshape(f_sigma,MX);
    
    % Antialiasing
    f_sigma(PX(1)+1,:)=0;
    f_sigma(:,PX(2)+1)=0;
    
    % % Influence of the complex brownian variance
    % f_sigma = 1/sqrt(prod(MX))*f_sigma;
    
    % % Choice to make the variance tensor explicitely independent of the
    % % resolution (i.e. independent of MX and of dX, but dependent on dX.*MX)
    % f_sigma = sqrt(prod(MX))*f_sigma;
    
    % Orthogonal gradient (from streamfunction to velocity)
    sigma_on_sq_dt(:,:,1)= 1i * ( - ky ) .* f_sigma;
    sigma_on_sq_dt(:,:,2)= 1i * ( + kx ) .* f_sigma;
else
    sigma_on_sq_dt = zeros([model.grid.MX 2]);
end


% Compute the bi-directional spectrum of sigma dBt
ft_sigma=abs(sigma_on_sq_dt).^2;
ft_sigma=sum(ft_sigma,3);

% % % Calcul of energy
% % % One has to divid by prod(model.grid.MX) because of the form of Parseval
% % % theorem for discrete Fourier transform
% % %trace_a_on_dt = sum(spectrum_a_sigma)
% % % trace_a_on_dt = 1/prod(model.grid.MX) * sum(spectrum_a_sigma);
% % % or equivalently
% % trace_a_on_dt = 1/prod(model.grid.MX) * sum(ft_sigma(:))
% trace_a = 1/prod(model.grid.MX) * sum(ft_sigma(:))

% Calcul of energy
% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
% One has to multiply by prod(model.grid.MX) because of the variance of
% the white in space noise
trace_a = 1/prod(model.grid.MX) * sum(ft_sigma(:));


if bool_plot
    
    % alpha = ( 3 - model.sigma.slope_absDif_sigma )/2;
    % eval(['norm_tr_a_theo = fct_norm_tr_a_theo_' ...
    %     model.sigma.type_spectrum '(model,k0,k_inf,alpha);']);
    % %norm_tr_a_theo = fct_norm_tr_a_theo(model,k0,k_inf,alpha);
    %
    % norm_tr_a_theo/trace_a_on_dt
    % % d_kappa
    % % prod(model.grid.MX)
    % % sqrt(prod(model.grid.MX))
    % % 1/(sqrt(prod(model.grid.MX))*d_kappa)
    % % d_kappa*prod(model.grid.MX)
    
    % Compute the omnidirectional spectrum of sigma dBt
    spectrum_a_sigma = idx' * ft_sigma(:);
    
    % Influence of the complex brownian variance
    spectrum_a_sigma = prod(MX)*spectrum_a_sigma;
    
    % % Choice to make the variance tensor explicitely independent of the
    % % resolution (i.e. independent of MX and of dX, but dependent on dX.*MX)
    % spectrum_a_sigma = (1/prod(MX))*spectrum_a_sigma;
    
    % Division by prod(model.grid.MX) because of the Parseval theorem for
    % discrete Fourier transform
    % Division by prod(model.grid.MX) again in order to the integration
    % of the spectrum over the wave number yields the energy of the
    % buoyancy averaged (not just integrated) over the space
    spectrum_a_sigma = 1/prod(model.grid.MX)^2 * spectrum_a_sigma;
    
    % % % Division by the wave number step
    % % spectrum_a_sigma = spectrum_a_sigma / d_kappa;
    %%
    
    % % % Division by prod(model.grid.MX) because of the Parseval theorem for
    % % % discrete Fourier transform
    % % Multiplication by prod(model.grid.MX) to normalize by the white-in-space noise
    % % variance
    % spectrum_a_sigma = prod(model.grid.MX) * spectrum_a_sigma;
    % %spectrum_a_sigma = 1/prod(model.grid.MX) * spectrum_a_sigma;
    
    % Division by the wave number step
    spectrum_a_sigma = spectrum_a_sigma / d_kappa;
    
    % % To remove the 2 pi which appear when we integrate k^(3-2-alpha) over the
    % % wave-vector angles and add the (2 pi)^2 which appear when we go from k to
    % % 2*pi*k
    % spectrum_a_sigma = spectrum_a_sigma * (2*pi);
    
    % Absolute diffusivity of the small-scale velocity (for test)
    trace_a_from_spectrum = d_kappa * sum(spectrum_a_sigma(:));
    
    
    %%
    
    % warning('debug here');
    % k2=(kx.^2+ky.^2);
    % reference_spectrum2 = k2 .^(1-alpha);
    % %reference_spectrum2 = k .^(2-2*alpha);
    % reference_spectrum2 = idx' * reference_spectrum2(:);
    % reference_spectrum2 = reference_spectrum2(2:end)';
    %
    % % Discretisation to go from continuous to discrete Fourier transform
    % reference_spectrum2 = 1/prod(model.grid.dX) * reference_spectrum2;
    %
    % % Division by prod(model.grid.MX) because of the Parseval theorem for
    % % discrete Fourier transform
    % reference_spectrum2 = 1/prod(model.grid.MX) * reference_spectrum2;
    %
    % % Division by the wave number step
    % reference_spectrum2 = reference_spectrum2 / d_kappa;
    %
    % % To remove the 2 pi which appear when we integrate k^(3-2-alpha) over the
    % % wave-vector angles
    % reference_spectrum2 = reference_spectrum2 * (2*pi);
    %
    % reference_spectrum(1)/reference_spectrum2(1)
    % reference_spectrum(1)/reference_spectrum2(1)*reference_spectrum2 ./ reference_spectrum
    %
    % figure(10);plot(kappa(2:end),reference_spectrum);hold on
    % plot(kappa(2:end),reference_spectrum2);
    
    %% Plot
    
    % % Make the plots appear at the same level thant the large-scale velocity spectrum
    % spectrum_a_sigma_plot = spectrum_a_sigma * mult_offset;
    spectrum_a_sigma_plot = spectrum_a_sigma;
    % reference_spectrum = reference_spectrum * mult_offset;
    
    taille_police = 12;
    
    X0 = [0 0];
    %     X0 = [10 20];
    
    %     figure10=figure(10);
    %     widthtemp = 12 ;
    %     heighttemp = 6;
    %     set(figure10,'Units','inches', ...
    %         'Position',[X0 widthtemp 2*heighttemp], ...
    %         'PaperPositionMode','auto');
    
    
    close(figure(10))
    figure10=figure(10);
    widthtemp = 12;
    heighttemp = 6;
    set(figure10,'Units','inches', ...
        'Position',[X0(1) X0(2) 2*widthtemp heighttemp], ...
        'PaperPositionMode','auto');
    
    
    %% Plot absolute diffusivity by scale
    
    subplot(1,2,2)
    
    loglog(kappa(2:end),reference_spectrum_a,'k')
    hold on;
    loglog(kappa(2:end),reference_spectrum_a_estim,'k--')
    name_plot = loglog(kappa(2:end),spectrum_w_a(2:end));
    set(name_plot,'LineWidth',LineWidth,...
        'MarkerSize',MarkerSize,...
        'Color',Color2);
    %     if nargin>2
    %         loglog(kappa(2:end),spectrum_w2_a(2:end),'c')
    %     end
    
    name_plot = loglog(kappa(2:end),spectrum_a_sigma_plot(2:end));
    set(name_plot,'LineWidth',LineWidth,...
        'MarkerSize',MarkerSize,...
        'Color',Color3);
    
    %     loglog(kappa(2:end),spectrum_w_a(2:end) + ...
    %         spectrum_a_sigma_plot(2:end),'g')
    ax=axis;
    ax(4)=max([spectrum_w_a ; ...
        reference_spectrum_a ; spectrum_a_sigma_plot]);
    % ax(4)=ax(4)*2;
    % ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w_a); ...
    %     max(spectrum_a_sigma_plot)]);
    d_ref = max([reference_spectrum_a(2:end) ; ...
        reference_spectrum_a_estim(2:end)]) ...
        / min([reference_spectrum_a(2:end) ;...
        reference_spectrum_a_estim(2:end)]);
    d_ref= sqrt(d_ref);
    mBound= min([reference_spectrum_a ; reference_spectrum_a_estim])...
        / d_ref;
    
    ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w_a); ...
        max(reference_spectrum_a_estim); max(reference_spectrum_a); ...
        max(spectrum_a_sigma_plot)]);
    ax(3) = min( [ax(3) ...
        mBound ...
        min(reference_spectrum_a) ...
        min(reference_spectrum_a_estim)]);
    %     ax(3) = min( [ax(3) ...
    %         min(spectrum_a_sigma_plot) ...
    %         min(reference_spectrum_a) ...
    %         min(reference_spectrum_a_estim)]);
    ax(1:2)=kappa([2 end]);
    if ax(4)>0
        axis(ax)
    end
    ax = axis;
    loglog(km_LS*[1 1],...
        [min(reference_spectrum_a_estim) ax(4)],'k--')
    loglog(...
        model.sigma.kappamin_on_kappamax_estim_slope * kappa(end)*[1 1],...
        [min(reference_spectrum_a_estim) ax(4)],'k--')
    loglog(model.sigma.kappamin_on_kappamax * kappa(end)*[1 1],...
        [min(reference_spectrum_a_estim) ax(4)],'k-.')
    hold off
    %     set(gca,'XGrid','on','XTickMode','manual');
    %     width = 9;
    %     height = 3;
    %     set(figure10,'Units','inches', ...
    %         'Position',[X0 2*width height], ...
    %         'PaperPositionMode','auto');
    %     set(gca,'YGrid','on')
    %     set(gca,...
    %         'Units','normalized',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    %     ylabel('$E(\kappa) \tau(\kappa) \bigl ( m^{3}.s^{-1}.{rad}^{-1} \bigr )$',...
    ylabel('$A(\kappa) = E(\kappa) \  \tau(\kappa) $',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    %     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    xlabel('$\kappa $',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'interpreter','latex',...
        'FontName','Times')
    title('ADSD',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    ax = axis;
    delta_ax = ax(4)-ax(3);
    
    %% Plot Spectrum
    subplot(1,2,1)
    idx_not_inf=~(isinf(log10(spectrum_w(2:end))) ...
        | spectrum_w(2:end)<1e-4*max(spectrum_w(2:end)) | isinf(kappa(2:end)));
    idx_not_inf = [false; idx_not_inf];
    line1= slope_ref * log10(kappa(2:end))  ;
    %     offset = -1 + mean(  log10(spectrum_w(idx_not_inf))  ...
    %         - line1(idx_not_inf(2:end)));
    offset = log10(spectrum_w(iii_k_LS(1))) ...
        - line1(iii_k_LS(1)-1);
    line1 = line1 + offset;
    ref=10.^line1;
    loglog(kappa(2:end),ref,'k');
    %loglog(kappa(2:end),ref,'--k');
    hold on;
    name_plot = loglog(kappa(2:end) , spectrum_w(2:end));
    set(name_plot,'LineWidth',LineWidth,...
        'MarkerSize',MarkerSize,...
        'Color',Color2);
    hold off;
    ax=axis;
    ax(4)=max([spectrum_w(2:end); ref(1:end)]);
    % ax(4)=ax(4)*2;
    
    min_ax= 10 ^(slope_ref * log10(kappa(2)*512/2) + offset) ;
    ax(3) = (model.odg_b/(1e-3))^2 * ...
        6e-2*(kappa(2)/kappa(end))*min([max(spectrum_w); max(ref)]);
    % ax(3)=6e-2*(kappa(2)/kappa(end))*min([max(spectrum_w); max(ref)']);
    ax(3) = min( [ax(3) min(ref) min_ax]);
    ax(1:2)=kappa(2)*[1 min(model.grid.MX)/2];
    if ax(4)>0
        axis(ax)
    end
    
    ax=axis;
    ax(3) = ax(4)-delta_ax;
    axis(ax);
    
    
    %     set(gca,'XGrid','on','XTickMode','manual');
    %     width = 4;
    %     height = 3;
    %     set(figure10,'Units','inches', ...
    %         'Position',[X0(1) X0(2) width height], ...
    %         'PaperPositionMode','auto');
    %     set(gca,'YGrid','on')
    
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    % ylabel('$|\hat{b}(\kappa)|^2$',...
    %     ylabel('$E(\kappa) \bigl ( m^{3}.s^{-2}.{rad}^{-1} \bigr )$',...
    ylabel('$E(\kappa) $',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    %     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    xlabel('$\kappa $',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'interpreter','latex',...
        'FontName','Times')
    title('Spectrum',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    %         'FontUnits','points',...
    %         'interpreter','latex',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    % %     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    %     xlabel('$\kappa $',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',taille_police,...
    %         'interpreter','latex',...
    %         'FontName','Times')
    %     title('Spectrum',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'interpreter','latex',...
    %         'FontSize',12,...
    %         'FontName','Times')
    
    
    
    drawnow
    
    %%
    
    subplot(1,2,2)
    set(gca,'XGrid','on','XTickMode','manual');
    %set(gca,'YGrid','on')
    drawnow
    subplot(1,2,1)
    set(gca,'XGrid','on','XTickMode','manual');
    %set(gca,'YGrid','on')
    drawnow
    
    %%
    
    subplot(1,2,2)
    width = 9;
    height = 3;
    set(figure10,'Units','inches', ...
        'Position',[X0 width height], ...
        'PaperPositionMode','auto');
    %set(gca,'YGrid','on')
    drawnow
    
    
    subplot(1,2,1)
    width = 9;
    height = 3;
    set(figure10,'Units','inches', ...
        'Position',[X0 width height], ...
        'PaperPositionMode','auto');
    %set(gca,'YGrid','on')
    
    
    subplot(1,2,2)
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    
    subplot(1,2,1)
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    
    %% Save plot
    drawnow
    
    folder_simu = model.folder.folder_simu;
    mkdir( [folder_simu ...
        '/AbsDiffByScale_sigma_dB_t_PostProcess']);
    eval( ['print -depsc ' folder_simu ...
        '/AbsDiffByScale_sigma_dB_t_PostProcess/'...
        day '.eps']);
    
    
%     folder_simu = [ pwd model.folder.folder_simu ];
%     folder_simu( folder_simu == '/' )='\';
%     mkdir( [folder_simu ...
%         '\AbsDiffByScale_sigma_dB_t_PostProcess']);
%     eval( ['print -depsc ' folder_simu ...
%         '\AbsDiffByScale_sigma_dB_t_PostProcess\'...
%         day '.eps']);
    
    % slope_w_a_comp_for_estim
    % km_LS
    % abs_diff_w
    % trace_a_on_2 = trace_a/2
end
end
function t = fct_unity_approx_(N_t)
% Approximation of unity
%

nx = 1:N_t;
nx=nx-mean(nx);
% see "New Numerical Results for the Surface Quasi-Geostrophic
% Equation", Constantin et al., J. Sci. Comput. (2012).
alpha = 36.;
%order = 30.;
order = 19.;
t = exp(-alpha*( (2./N_t).*abs(nx) ).^order);


% sslop=8;
% t=ones(1,N_t);
% t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

end
% function t = fct_unity_approx_old(N_t)
% % Approximation of unity
% %
%
% sslop=8;
% t=ones(1,N_t);
% t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
%
% end