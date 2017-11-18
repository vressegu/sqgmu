function [slope,km]=fct_estim_spectrum_slope(model,f)
% This function creates some plot online and save it
%

bool_plot = false;

%% Get paramters
if isfield(model,'spectrum_theo')
    model=rmfield(model,'spectrum_theo');
end
iii = isinf(abs(f));
f = f - mean(f(~iii));

% Grid
x = model.grid.x;
y = model.grid.y;
[X,Y]=ndgrid(x,y);

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
% folder_simu = model.folder.folder_simu;
plot_moments = false;
map = model.folder.colormap;

%% One particle
X0=[0 0];

% if model.mirror
%     f=f(:,1:model.grid.MX(2)/2,:,:);
%     y=y(1:model.grid.MX(2)/2);
% end

width = 3.3;

%     height = 3.2;

ax = [x(end)-x(1) y(end)-y(1)] ;
aspect_ratio = ax(2)/ax(1);
height = aspect_ratio * width;

%% Cleaning

lon = model.grid.lonlat.lon;
lat = model.grid.lonlat.lat;
[LON,LAT]=ndgrid(lon,lat);

if bool_plot
    figure;imagesc(x,y,f');axis xy;axis equal;
end

f_sph = fct_cart2sph(x,y,f,LON,LAT,model.grid.lonlat.lonlat_ref);

if bool_plot
    figure;imagesc(lon,lat,f_sph');axis xy;axis equal;
end

% s = size(f_sph);
% f_sph(iii)=0;
% f_sph = reshape(f_sph,s);

marge_lon =0;
marge_lat =0;
crop_f = f_sph;
while ...
        marge_lon < ceil(min(length(lon))/2) ...
        && ...
        marge_lat < ceil(min(length(lat))/2) ...
        && ...
        any(isnan(crop_f(:)))
    marge_lon = marge_lon+1;
    crop_f = crop_f(2:end-1,:);
    % To cropped latitutes slower
    if mod(marge_lon,2)==0
        marge_lat = marge_lat+1;
        crop_f = crop_f(:,2:end-1);
    end
    if bool_plot
        marge_lon
        marge_lat
        figure(3);imagesc(crop_f);
        drawnow;pause(0.01);
    end
end
remove = 10;
slop_size_ratio=3;
sslop=ceil(size(LON)/slop_size_ratio);
marge_lon = max([0 ceil(marge_lon-sslop(1)/remove)]);
marge_lat = max([0 ceil(marge_lat-sslop(2)/remove)]);
% marge_lon = max([0 ceil(marge_lon-sslop(1)/4)]);
% marge_lat = max([0 ceil(marge_lat-sslop(2)/4)]);
% marge_lon = ceil(marge_lon-sslop(1)/4);
% marge_lat = ceil(marge_lat-sslop(2)/4);

% warning('bidouille')
% marge(2)=marge(2)/2;

% iii = isnan(f_sph) | isinf(abs(f_sph));
% iiilon=find(~any(iii,2));
% id_iiilon = iiilon([1 end]);
% len_iiilon = 1 + id_iiilon(2)-id_iiilon(1);
% mask_lon = [ zeros(1,id_iiilon(1)-1) ...
%              fct_unity_approx6(len_iiilon) ...
%              zeros(1,length(model.grid.lonlat.lon)-id_iiilon(2)) ];
% iiilat=find(~any(iii,1));
% id_iiilat = iiilat([1 end]);
% len_iiilat = 1 + id_iiilat(2)-id_iiilat(1);
% mask_lat = [ zeros(1,id_iiilat(1)-1) ...
%              fct_unity_approx6(len_iiilat) ...
%              zeros(1,length(model.grid.lonlat.lat)-id_iiilat(2)) ];

% marge
% mask_lon = [ zeros(1,marge) ...
%              fct_unity_approx6(...
%              length(model.grid.lonlat.lon)-2*marge) ...
%              zeros(1,marge) ];
% mask_lat = [ zeros(1,marge) ...
%              fct_unity_approx6(...
%              length(model.grid.lonlat.lat)-2*marge) ...
%              zeros(1,marge) ];
% mask_lon = [ zeros(1,max([0 marge(1)])) ...
%              fct_unity_approx6(...
%              length(model.grid.lonlat.lon)-2*max([0 marge(1)])) ...
%              zeros(1,max([0 marge(1)])) ];
% mask_lat = [ zeros(1,max([0 marge(2)])) ...
%              fct_unity_approx6(...
%              length(model.grid.lonlat.lat)-2*max([0 marge(2)])) ...
%              zeros(1,max([0 marge(2)])) ];
mask_lon = [ zeros(1,marge_lon) ...
    fct_unity_approx6(...
    length(model.grid.lonlat.lon)-2*marge_lon) ...
    zeros(1,marge_lon) ];
mask_lat = [ zeros(1,marge_lat) ...
    fct_unity_approx6(...
    length(model.grid.lonlat.lat)-2*marge_lat) ...
    zeros(1,marge_lat) ];

% Remove boundaries
mask_boundaries = mask_lon' * mask_lat;

% Remove boundaries
% mask_boundaries = ...
%     fct_unity_approx6(length(model.grid.lonlat.lon))' * ...
%     fct_unity_approx6(length(model.grid.lonlat.lat));
f_sph = bsxfun(@times, mask_boundaries, f_sph);

if bool_plot
    figure;imagesc(lon,lat,f_sph');axis xy;axis equal;
end

f_cart = fct_sph2cart(lon,lat,f_sph,X,Y,model.grid.lonlat.lonlat_ref);
if bool_plot
    figure;imagesc(x,y,f_cart');axis xy;axis equal;
end

iii = isnan(f_cart) | isinf(abs(f_cart));
s = size(f_cart);
f_cart(iii)=0;
f_cart = reshape(f_cart,s);

if bool_plot
    figure;imagesc(x,y,f_cart');axis xy;axis equal;
end

%% Remove zero padding
iii = (f_cart~=0);
iiilon=find(any(iii,2));
id_iiilon = iiilon([1 end]);
len_iiilon = 1 + id_iiilon(2)-id_iiilon(1);
if mod(len_iiilon,2)==1
    id_iiilon(2)=id_iiilon(2)-1;
end
% mask_lon = [ zeros(1,id_iiilon(1)-1) ...
%              fct_unity_approx6(len_iiilon) ...
%              zeros(1,length(model.grid.lonlat.lon)-id_iiilon(2)) ];
iiilat=find(any(iii,1));
id_iiilat = iiilat([1 end]);
len_iiilat = 1 + id_iiilat(2)-id_iiilat(1);
if mod(len_iiilat,2)==1
    id_iiilat(2)=id_iiilat(2)-1;
end
% mask_lat = [ zeros(1,id_iiilat(1)-1) ...
%              fct_unity_approx6(len_iiilat) ...
%              zeros(1,length(model.grid.lonlat.lat)-id_iiilat(2)) ];
f_cart = f_cart(id_iiilon(1):id_iiilon(2),id_iiilat(1):id_iiilat(2));
x = x(id_iiilon(1):id_iiilon(2));
y = y(id_iiilat(1):id_iiilat(2));
model.gridref = model.grid;
model.grid.MX = size(f_cart);

if bool_plot
    figure(30);
    imagesc(x,y,f_cart');axis xy;axis equal;
    title('plot for spectrum')
end

% % Remove boundaries
% mask_boundaries = ...
%     fct_unity_approx6(length(model.grid.lonlat.lon))' * ...
%     fct_unity_approx6(length(model.grid.lonlat.lat));
% % mask_boundaries = ...
% %     fct_unity_approx6(model.grid.MX(1))' * ...
% %     fct_unity_approx6(model.grid.MX(2));
%
% mask_boundaries = ( ...
%     fct_sph2cart(model.grid.lonlat.lon,model.grid.lonlat.lon,...
%     mask_boundaries,X,Y,model.grid.lonlat.lonlat_ref, ...
%     'linear') == 1 );
% % % mask_keep_cart = fct_sph2cart(lon_global,lat_global,mask_keep_,X,Y,lonlat_ref);
% % s = size(mask_keep_cart);
%
%
% f = bsxfun(@times, mask_boundaries, f);

fft_w = fft2(f_cart);

%% PLots

if ~(isfield(model.folder,'hide_plots'))
    model.folder.hide_plots = false;
end
if ~model.folder.hide_plots;
    
    X0=[3.3 1];
    % close(figure(4))
    figure4=figure(4);
    % figure4=figure;
    
    widthtemp = 12;
    heighttemp = 6;
    set(figure4,'Units','inches', ...
        'Position',[X0(1) X0(2) widthtemp heighttemp], ...
        'PaperPositionMode','auto');
end
[slope,km]=fct_spectrum_estim_slope( model,fft_w(:,:,:,id_part),'b');
% fct_spectrum( model,fft_w(:,:,:,id_part),'b');
if ~model.folder.hide_plots;
    set(gca,'XGrid','on','XTickMode','manual');
    width = 4.5;
    % width = 4;
    height = 3;
    set(figure4,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    set(gca,'YGrid','on')
    
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('$\overline{\Gamma}_T$',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
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
    drawnow
    % eval( ['print -depsc ' folder_simu '/Spectrum/' day '.eps']);
    
end

end

