function [slope,km,spectrum,name_plot] = fct_spectrum_estim_slope(model,ft,color)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end

% % figure;plot(abs((ft(1,:))));
% figure;plot(abs([ft(1,:) ft(1,:)]));
% % figure;plot(real(ifft(ft(1,:))));
% % figure;imagesc(real(ifft2(ft')));axis equal;axis xy

% Square modulus of the Fourier transform
ft=abs(ft).^2;

% Get parameters
MX=model.grid.MX;
PX=MX/2;
dX=model.grid.dX;
if any(size(ft)~=MX)
    error('wrong size');
end
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

% ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% ky=2*pi/model.grid.dX(2)*ky;
% % figure;imagesc(ft');axis xy
% fty=ft(1,:);
% % fty = [fty fty];
% figure;plot(ky,(fty),'o')
% hold on;
% Ly = model.grid.dX(2)*model.grid.MX(2);
% warning('specification of sigma again');
% sigma = 0.002* Ly;
% % sigma = 0.01 * Ly;
% % % sigma = 0.03 * Ly;
% fty_theo = 1/prod(model.grid.dX) ...
%     * model.grid.MX(2) * Ly * (2*pi)^3 ...
%     * sigma^2 * exp( - sigma^2 * ky.^2);
% s = sort(fty);
% s_theo = sort(fty_theo);
% s = s(end-1);
% s_theo = s_theo(end-1);
% rate = s/s_theo;
% fty_theo = rate * fty_theo;
% plot(ky,(fty_theo),'r.')
% hold off
% % figure;plot(fty)

persistent idxref MXref dXref kidx

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    
    % Remove aliasing
    ft(PX(1),:)=0;
    ft(:,PX(2))=0;
    
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

%% Spectrum
% Integration over the rings of iso wave number
spectrum = idxref' * ft(:);

% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/prod(model.grid.MX)^2 * spectrum;

% Division by the wave number step
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

%% Theoritical spectrum tail
if isfield(model,'spectrum_theo')
    %     spectrum_theo = model.spectrum_theo.coef1 * kidx ...
    spectrum_theo = model.spectrum_theo.coef1 ...
        .* exp( -1/2 * model.spectrum_theo.coef2 * kidx .^2 );
    
    % Division by prod(model.grid.MX) in order to the integration
    % of the spectrum over the wave number yields the energy of the
    % buoyancy averaged (not just integrated) over the space
    %     spectrum_theo = 1/prod(model.grid.MX) * spectrum_theo;
    spectrum_theo = 1/prod(model.grid.dX.*model.grid.MX) * spectrum_theo;
    
    spectrum_theo = 1/(2*pi)^2 * spectrum_theo;
    
    %     spectrum_theo = 1/prod(model.grid.MX) * spectrum_theo;
    %     spectrum_theo = 1/prod(model.grid.MX) * spectrum_theo;
end

% spectrum_theo = spectrum_theo ./ kidx;
% spectrum = spectrum ./ kidx';

%% Plot
if isfield(model,'spectrum_theo') && ~model.folder.hide_plots
    loglog(kidx(2:end),spectrum_theo(2:end),'r');
    hold on;
end
if isfield(model,'gridref') && any(model.gridref.MX > model.grid.MX)
    M_kappa=min(model.gridref.MX);
    km = 2*pi*max(1./model.gridref.dX)/M_kappa;
else
    km = kidx(2);
end
% loglog(kidx(2:end),ref,'--k');
% hold on;
if ~model.folder.hide_plots
    name_plot = loglog(kidx(2:end) , spectrum(2:end) ,color);
    ax=axis;
    ax(4)=max(spectrum(2:end));
    % ax(4)=max([spectrum(2:end); ref']);
    % min_ax= 10 ^(-5/3 * log10(kidx(2)*512/2) + offset) ;
    ax(3)=min(spectrum(2:end));
    % ax(3)=6e-2*(kidx(2)/kidx(end))*min([max(spectrum); max(ref)']);
    % ax(3) = min( [ax(3) min(ref) min_ax]);
    % ax(3) = max( [ax(3) ax(4)*1e-6 ]);
    % % ax(3) = max( [ax(3) ax(4)*1e-10 ]);
    % ax(1:2)=kidx(2)*[1 min(model.grid.MX)/2];
    % if ax(4)>0
    %     axis(ax)
    % end
    ax(3) = max( [ax(3) ax(4)*1e-7 ]);
    % ax(3) = max( [ax(3) ax(4)*1e-10 ]);
    ax(1:2)=[km kidx(2)*min(model.grid.MX)/2];
    if ax(4)>0
        axis(ax)
    end
end

%% Estimation of slope
[~,i_first]=max(spectrum(2:end));
i_first = i_first +1;
iii = i_first:length(spectrum);
spectrum = spectrum(iii);
kkk=kidx(iii);
km=kkk(1);
iii=~(isinf(abs(spectrum))|isnan(spectrum)|(spectrum<=eps)|...
    (kkk'>kidx(end)/6));
%     (kkk'>kidx(end)/3));
%     (kkk'>10*km));
logspectrum = log10(spectrum(iii));
logkkk=log10(kkk(iii));
Z = logkkk';
Z(:,2) = 1;
beta = Z\logspectrum;
slope = beta(1);
offset = beta(2);

% slope = 1/(logkkk * logkkk') * logspectrum * logkkk';

%% Plot estim
if ~model.folder.hide_plots
    % idx_not_inf=~(isinf(log10(spectrum(2:end))) ...
    %     | spectrum(2:end)<1e-4*max(spectrum(2:end)) | isinf(kidx(2:end)'));
    % line1= slope * log10(kidx(2:end))  ;
    line1= slope * log10(kidx(2:end)) + offset ;
    % offset = -1 + mean(  log10(spectrum(idx_not_inf)')  - line1(idx_not_inf));
    % line1 = line1 + offset;
    ref=10.^line1;
    hold on;
    loglog(kidx(2:end),ref,'--k');
    hold off
end
