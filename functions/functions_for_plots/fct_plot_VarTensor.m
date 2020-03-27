function fct_plot_VarTensor(model,day)
% This function creates some plot online and save it
%

%% Get paramters

% Grid
x = model.grid.x;
y = model.grid.y;
My = model.grid.MX(2);

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
plot_moments = model.advection.plot_moments;
%plot_epsilon_k = model.advection.plot_epsilon_k;
map = model.folder.colormap;

%% One particle
X0=[0 0];
VarTensor = 2 * model.advection.coef_diff(:,:,:,id_part);
VarTensor = bsxfun( @times, ones(model.grid.MX), VarTensor) ;
VarTensor = sum(VarTensor,3);
% T_adv_part = real(ifft2( fft_b_adv_part(:,:,1,id_part) ));

% if ( (eval(day) == 0) && ...
%         strcmp(model.type_data,'Perturbed_vortices') )
%     width = 3.2;
%     height = 3.2;
%     figure1=figure(1);
%     set(figure1,'Units','inches', ...
%         'Position',[X0(1) X0(2) width height], ...
%         'PaperPositionMode','auto');
%     contourf(x,y,T_adv_part');
%     x= model.grid.dX(1)*(0:model.grid.MX(1)-1);
%     y= model.grid.dX(2)*(0:model.grid.MX(2)-1);
%     Lx = model.grid.dX(1) * model.grid.MX(1);
%     sigma= 2 * Lx/15;
%     center1x=x(1/4*model.grid.MX(1)+1);
%     center1y=y(1/4*model.grid.MX(2)+1);
%     nsig=40;
%     dist = 1.5;
%     rate = 0.3;
%     sigma = sigma/nsig;
%     center1x= 2e4;
%     center1y=y(1/4*model.grid.MX(2)+1);
%     coord1=[center1x center1y];
%     size_square = 10e4;
%     redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
%     hold on;
%     if strcmp(model.type_data,'Perturbed_vortices')
%         plot(redline1(:,1),redline1(:,2),'r','LineWidth',3);
%     elseif strcmp(model.type_data,'spot6')
%         plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
%     else
%         error('wrong type of data?');
%     end
%     center1x= x(1/2*model.grid.MX(1)+1) - 2e4;
%     center1y=y(3/4*model.grid.MX(2)+1);
%     coord1=[center1x center1y];
%     redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
%     plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
%     hold off;
%
% else
width = 3.3;
height = 3.2;
figure772=figure(772);
close(figure772);
figure772=figure(772);
set(figure772,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
imagesc(x,y,VarTensor');
% end
% caxis([-1 1]*model.odg_b);
% % caxis([-1 1]*1e-3);
cax = caxis; cax(1)=0;caxis(cax);
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
% % title({'One realization', ...
% %     ['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
% title(['\hspace{0.5cm} $t=' num2str(day) '$ day '],...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times')
axis xy; axis equal
% colormap(map)
colorbar
drawnow
mkdir([folder_simu '/VarTensor']);
eval( ['print -depsc ' folder_simu '/VarTensor/'...
    num2str(day) '.eps']);


