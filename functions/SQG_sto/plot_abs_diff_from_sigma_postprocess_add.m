function plot_abs_diff_from_sigma_postprocess_add(model,...
    fft_sigma_dBt_on_sq_dt,color)
% - fft_sigma_dBt_on_sq_dt is the Fourier transform of 
% the kernel \tilde sigma up to a multiplicative constant
%

%%


% id_part = 1;
% 
% fft_sigma_dBt_on_sq_dt = fft_sigma_dBt_on_sq_dt(:,:,:,id_part);


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
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;
fft_sigma_dBt_on_sq_dt(PX(1)+1,:,:,:)=0;
fft_sigma_dBt_on_sq_dt(:,PX(2)+1,:,:)=0;

% Compute the bi-directional spectrum of sigma dBt
ft_sigma=abs(fft_sigma_dBt_on_sq_dt).^2;
ft_sigma=sum(ft_sigma,3);
ft_sigma=mean(ft_sigma,4);

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
idx = sparse( bsxfun(@le, kappa, k ) );
idx = idx & sparse( bsxfun(@lt,k, [ kappa(2:end) kappa(end)+d_kappa ] ) );

kappa = kappa';

%% Spectrum


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
trace_a = 1/prod(model.grid.MX)^2 * sum(ft_sigma(:));



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

% % Influence of the complex brownian variance
% spectrum_a_sigma = prod(MX)*spectrum_a_sigma;

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
trace_a_from_spectrum = d_kappa * sum(spectrum_a_sigma(:))

%%

% Absolute diffusivity by scale with theoretical slope
reference_spectrum_a = kappa(2:end) .^ slope_ref_a ;
%reference_spectrum_a = kappa(2:end) .^ model.sigma.slope_sigma ;

iii_k_LS = ceil(1/4*length(spectrum_a_sigma ));

% Setting offset
mult_offset_spectrum_a = spectrum_a_sigma (iii_k_LS(1)) ...
    / reference_spectrum_a(iii_k_LS(1)-1) ;
reference_spectrum_a = mult_offset_spectrum_a * reference_spectrum_a;

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
X0 = [10 20];

%     figure10=figure(10);
%     widthtemp = 12 ;
%     heighttemp = 6;
%     set(figure10,'Units','inches', ...
%         'Position',[X0 widthtemp 2*heighttemp], ...
%         'PaperPositionMode','auto');


% close(figure(10))
figure10=figure(10);
subplot(1,2,2)
widthtemp = 12;
heighttemp = 6;
set(figure10,'Units','inches', ...
    'Position',[X0(1) X0(2) 2*widthtemp heighttemp], ...
    'PaperPositionMode','auto');

%% Plot absolute diffusivity by scale

hold on;
% loglog(kappa(2:end),reference_spectrum_a,'k')

% loglog(kappa(2:end),reference_spectrum_a_estim,'k--')
%     name_plot = loglog(kappa(2:end),spectrum_a_sigma (2:end));
%     set(name_plot,'LineWidth',LineWidth,...
%         'MarkerSize',MarkerSize,...
%         'Color',Color2);
%     %     if nargin>2
%     %         loglog(kappa(2:end),spectrum_w2_a(2:end),'c')
%     %     end

name_plot = loglog(kappa(2:end),spectrum_a_sigma_plot(2:end));
set(name_plot,'LineWidth',LineWidth,...
    'MarkerSize',MarkerSize,...
    'Color',color);

%     loglog(kappa(2:end),spectrum_a_sigma (2:end) + ...
%         spectrum_a_sigma_plot(2:end),'g')
ax=axis;
ax(4)=max([ax(4) ; spectrum_a_sigma  ; ...
    reference_spectrum_a ; spectrum_a_sigma_plot]);
% % ax(4)=ax(4)*2;
% % ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_a_sigma ); ...
% %     max(spectrum_a_sigma_plot)]);
% d_ref = max([reference_spectrum_a(2:end) ; ...
%     reference_spectrum_a_estim(2:end)]) ...
%     / min([reference_spectrum_a(2:end) ;...
%     reference_spectrum_a_estim(2:end)]);
% d_ref= sqrt(d_ref);
% mBound= min([reference_spectrum_a ; reference_spectrum_a_estim])...
%     / d_ref;

ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_a_sigma ); ...
    max(reference_spectrum_a); ...
    max(spectrum_a_sigma_plot)]);
ax(3) = min( [ax(3) ...
    min(reference_spectrum_a) ]);
% ax(3) = min( [ax(3) ...
%     mBound ...
%     min(reference_spectrum_a) ...
%     min(reference_spectrum_a_estim)]);
% %     ax(3) = min( [ax(3) ...
% %         min(spectrum_a_sigma_plot) ...
% %         min(reference_spectrum_a) ...
% %         min(reference_spectrum_a_estim)]);
ax(1:2)=kappa([2 end]);
if ax(4)>0
    axis(ax)
end
ax = axis;
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


% %     set(gca,'XGrid','on','XTickMode','manual');
% %     width = 4;
% %     height = 3;
% %     set(figure10,'Units','inches', ...
% %         'Position',[X0(1) X0(2) width height], ...
% %         'PaperPositionMode','auto');
% %     set(gca,'YGrid','on')
% 
% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',taille_police,...
%     'FontName','Times')
% % ylabel('$|\hat{b}(\kappa)|^2$',...
% %     ylabel('$E(\kappa) \bigl ( m^{3}.s^{-2}.{rad}^{-1} \bigr )$',...
% ylabel('$A(\kappa) $',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',taille_police,...
%     'FontName','Times')
% %     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
% xlabel('$\kappa $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',taille_police,...
%     'interpreter','latex',...
%     'FontName','Times')
% title('ADSD',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times')
% %         'FontUnits','points',...
% %         'interpreter','latex',...
% %         'FontSize',taille_police,...
% %         'FontName','Times')
% % %     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
% %     xlabel('$\kappa $',...
% %         'FontUnits','points',...
% %         'FontWeight','normal',...
% %         'FontSize',taille_police,...
% %         'interpreter','latex',...
% %         'FontName','Times')
% %     title('Spectrum',...
% %         'FontUnits','points',...
% %         'FontWeight','normal',...
% %         'interpreter','latex',...
% %         'FontSize',12,...
% %         'FontName','Times')



drawnow

%%

%subplot(1,2,2)
set(gca,'XGrid','on','XTickMode','manual');
set(gca,'YGrid','on')
%%

width = 9;
height = 3;
set(figure10,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
%set(gca,'YGrid','on')

set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
drawnow

%% Save plot
drawnow

folder_simu = model.folder.folder_simu;
eval( ['print -depsc ' folder_simu ...
    '/AbsDiffByScale_sigma_dB_t_PostProcess_PostProcess.eps']);

% slope_w_a_comp_for_estim
% km_LS
% abs_diff_w
% trace_a_on_2 = trace_a/2

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