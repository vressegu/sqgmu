function [spectrum,name_plot] = fct_spectrum53_HR_error(model,ft,color)
% Compute the spectrum of a function
%

if nargin < 3
    color='b';
end

ft=abs(ft).^2;

MX=model.grid.MX;
dX=model.grid.dX;

if any(size(ft)~=MX)
    error('wrong size');
end
% [Mx,My]=size(ft);
% MX=[Mx My];
% M=floor(sqrt(Mx*My));
% M=min(MX);
% M_kappa=min(MX/4);
% M=min(MX/8);
persistent idxref MXref dXref kidx

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    %% Define wave number
    if any( mod(MX,2)~=0)
        error('the number of grid points by axis need to be even');
    end
    PX=MX/2;
    ft(PX(1),:)=0;
    ft(:,PX(2))=0;
    
    
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    kx=2*pi/model.grid.dX(1)*kx;
    ky=2*pi/model.grid.dX(2)*ky;
    [kx,ky]=ndgrid(kx,ky);
    k=sqrt(kx.^2+ky.^2);
    k(PX(1)+1,:)=inf;
    k(:,PX(2)+1)=inf;
    k=k(:);
    
%     kx=1/(MX(1)*dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
%     ky=1/(MX(2)*dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% %     kx=1/Mx*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
% %     ky=1/My*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
%     [kx,ky]=ndgrid(kx,ky);
% %     k=log10(kx.^2+ky.^2) ; clear kx ky
%     k=log10(kx.^2+ky.^2) /2; clear kx ky
%     k(PX(1),:)=inf;
%     k(:,PX(2))=inf;
%     % /2 is for the sqrt of the norm 2
% %     % k=sqrt(kx.^2+ky.^2); clear kx ky
%     k=k(:);
%     % k=repmat(k,[1 M]);
    
    %% Order it
    
    M_kappa=min(MX);
    P_kappa= M_kappa/2;
    d_kappa = max(1./dX);
    kidx=1/(M_kappa)* (0:(P_kappa-1)) ;
    kidx=2*pi*d_kappa*kidx;
    %% Masks associated with the rings of iso wave number
    d_kappa = kidx(2) - kidx(1);
    idx = sparse( bsxfun(@le,kidx, k ) );
    idx = idx & sparse( bsxfun(@lt,k, [ kidx(2:end) kidx(end)+d_kappa ] ) );
    idxref=idx;
end

% for j=1:P_kappa-1
%     imagesc(reshape(idxref(:,j),MX))
%     drawnow;
%     pause
% end

%% Spectrum
% ft = 1/prod(model.grid.MX)^2 * ft;
spectrum = idxref' * ft(:);
% spectrum = [ft(1); spectrum]; % Add the constant field

% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration 
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/prod(model.grid.MX)^2 * spectrum;

% Division by the wave number step
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

%% Plot
idx_not_inf=~(isinf(log10(spectrum(2:end)))| ...
    spectrum(2:end)<1e-4*max(spectrum(2:end)) | isinf(kidx(2:end)'));
% idx_not_inf=~(isinf(log10(spectrum))| spectrum<1e-4*max(spectrum) | isinf(kidx'));
line1= -5/3 * log10(kidx)  ;
idx_not_inf=[false ; idx_not_inf];
% offset = -8.5 - (-5/3) * (-5)  ;
offset = -1 + mean(  log10(spectrum(idx_not_inf)')  - line1(idx_not_inf));
% offset = 1+log10(3) - (-5/3) * (-5) - 2*log10(prod(model.grid.MX)/(256^2)) ;
% % offset = mean(  log10(spectrum(idx_not_inf)')  - line1(idx_not_inf));
line1 = line1 + offset;
ref=10.^line1(2:end);
% loglog(kidx(2:end),ref,'--k');
% hold on;
name_plot = loglog(kidx(2:end) , spectrum(2:end) ,color);
hold on;
% loglog(kidx(2:end) , spectrum(2:end) )
ax=axis;
% ax(4)=max([spectrum; ]);
% ax(3)=6e-2*(kidx(2)/kidx(end))*min([max(spectrum);]);

ax(1:2)=kidx(2)*[1 128/2];
% ax(1:2)=kidx(2)*[1 256/2];
% % ax(1:2)=kidx(2)*[1 512/2];
% % % ax(1:2)=kidx([2 end]);
if max(spectrum(2:end))>ax(4)
    ax(4)=max(spectrum(2:end));
end
if ax(4)>0
    axis(ax)
end
