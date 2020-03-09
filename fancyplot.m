% This script aims to produce journal-quality plots of each step in the Solar PV project.

% Startup
clc
u = symunit;
buildings = ["ADMIN"; "BioSci";  "CAM"; "CCIS"; "DICE"; ...
             "ETLC";     "GSB"; "NINT"; "NREF";  "SUB"];
         
% for i = 1:10
    i = 4;
    length = 1600;
    folder = "output/fancyplot/";

    name = buildings(i);
    b = Building(name,length);

    % Save for Nima
%     eval(name+"_table = b.info");
%     eval(name+"_pix2m = b.pix2m");
%     save("output/transfer/"+name+"_info.mat",name+"_table",name+"_pix2m");
%     eval("table = b.info");
%     eval("pix2m = b.pix2m");
%     save("output/transfer/"+name+"_info.mat","table","pix2m");
% end

%% ~~~~~~~~~~~~~~
% ~~ Plotting ~~
% ~~~~~~~~~~~~~~

fontsize = 14;

% Height map
figure(1);
b.plot(b.Z);
ax = gca; ax.FontSize = fontsize;
saveas(figure(1),folder+"maps/"+name+"_hmap.png");

% Height+Terrain map
hold on; alpha 0.5; surf(b.X,b.Y,b.terrain); hold off;
ax = gca; ax.FontSize = fontsize;
view(3);
saveas(figure(1),folder+"maps/"+name+"_hmap_terrain.png");
close(1);

% Colour map
figure(1);
b.plot(b.map_r);
ax = gca; ax.FontSize = fontsize;
saveas(figure(1),folder+"maps/"+name+"_cmap.png");
close(1);

% Height+Colour maps combined
figure(1);
width=1920;height=1080;
set(gcf,'position',[0,0,width,height]);
subplot(1,2,1);
b.plot(b.Z);
ax = gca; ax.FontSize = fontsize;
subplot(1,2,2);
b.plot(b.map_r);
ax = gca; ax.FontSize = fontsize;
saveas(figure(1),folder+"maps/"+name+"_both.png");
close(1);

% Surfaces
figure(1);
b.plot(b.Z .* b.issurface);
ax = gca; ax.FontSize = fontsize;
title("Surfaces of "+name);
saveas(figure(1),folder+"identified/"+name+"_surfaces.png");
close(1);

% Roofs
figure(1);
b.plot(b.Z .* b.isroof);
ax = gca; ax.FontSize = fontsize;
title("Rooftops of "+name);
saveas(figure(1),folder+"identified/"+name+"_roofs.png");
close(1);

% Surface/Roof Comparison
figure(1);
% ax1 = subplot(1,3,1);
% b.plot(b.Z .* b.issurface + 0.10); colorbar off; title("Surfaces"+newline+"(not edge, ground, or tiny)");
% ax1.FontSize = fontsize; 
% hold on;
% surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
% hold off;
% view(3); camlight('headlight'); view(2);
ax2 = subplot(1,2,1);
% b.plot(b.Z .* b.isroof + 0.10); colorbar off; title("Roofs"+newline+"(surfaces that can fit >1 SPV)");
b.plot(b.Z .* b.isroof + 0.10); colorbar off; title("Roofs"+newline+"(not edge, ground, or tiny)");
hold on;
surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
ax2.FontSize = fontsize;
hold off;
view(3); camlight('headlight'); view(2);
ax3 = subplot(1,2,2);
b.plot(b.Z .* b.isroof_viable + 0.10); colorbar off; title("Viable roofs"+newline+"(roofs that are flat/slanted)");
hold on;
surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
ax3.FontSize = fontsize;
hold off;
view(3); camlight('headlight'); view(2);
width=1920;height=1080;
set(gcf,'position',[0,0,width,height]);
saveas(figure(1),folder+"comparisons/"+name+"_identified.png");
close(1);

% Roof types, labeled
figure(1);
rooftype = zeros(size(b.roofs)); % Will be labeled according to type
for i = 1:max(b.roofs,[],'all')
    type = b.info.type(b.info.roof==i);
    rooftype(b.roofs==i) = 3*(type=='flat') + 2*(type=='slanted') + 1*(type=='irregular');
end
imagesc(rooftype); axis image;
% cmap = cat(1,0.15*[1,1,1],pink(3));
cmap = cat(1,0.15*[1,1,1],[0.5774,0,0],[0.8165    0.7165    0.3774],[0.9,0.9,0.8]);
colormap(cmap);
L = line(ones(size(cmap,1)-1),ones(size(cmap,1)-1),'LineWidth',10);
set(L,{'color'},mat2cell(cmap(2:end,:),ones(1,size(cmap,1)-1),size(cmap,2)));            % set the colors according to cmap
legend(L, {'irregular', 'slanted','flat'}, 'Location','southwest')
ax = gca; ax.FontSize = fontsize*2/3;
saveas(figure(1),folder+"maps/"+name+"_roof_types.png");
close(1);

%% LOOP THROUGH ALL BUILDINGS

for i = 1:size(buildings,1)
    name = buildings(i)
    b = Building(name,length);
    eval("table = b.info");
    eval("pix2m = b.pix2m");
    save("output/transfer/"+name+"_info.mat","table","pix2m");
    
%     % Plotting
%     % Height map
%     figure(1);
%     b.plot(b.Z);
%     ax = gca; ax.FontSize = fontsize;
%     saveas(figure(1),folder+"maps/"+name+"_hmap.png");
% 
%     % Height+Terrain map
%     hold on; alpha 0.5; surf(b.X,b.Y,b.terrain); hold off;
%     ax = gca; ax.FontSize = fontsize;
%     view(3);
%     saveas(figure(1),folder+"maps/"+name+"_hmap_terrain.png");
%     close(1);
% 
%     % Colour map
%     figure(1);
%     b.plot(b.map_r);
%     ax = gca; ax.FontSize = fontsize;
%     saveas(figure(1),folder+"maps/"+name+"_cmap.png");
%     close(1);
% 
%     % Height+Colour maps combined
%     figure(1);
%     width=1920;height=1080;
%     set(gcf,'position',[0,0,width,height]);
%     subplot(1,2,1);
%     b.plot(b.Z);
%     ax = gca; ax.FontSize = fontsize;
%     subplot(1,2,2);
%     b.plot(b.map_r);
%     ax = gca; ax.FontSize = fontsize;
%     saveas(figure(1),folder+"maps/"+name+"_both.png");
%     close(1);
% 
%     % Surfaces
%     figure(1);
%     b.plot(b.Z .* b.issurface);
%     ax = gca; ax.FontSize = fontsize;
%     title("Surfaces of "+name);
%     saveas(figure(1),folder+"identified/"+name+"_surfaces.png");
%     close(1);
% 
%     % Roofs
%     figure(1);
%     b.plot(b.Z .* b.isroof);
%     ax = gca; ax.FontSize = fontsize;
%     title("Rooftops of "+name);
%     saveas(figure(1),folder+"identified/"+name+"_roofs.png");
%     close(1);
% 
%     % Surface/Roof Comparison
%     figure(1);
%     % ax1 = subplot(1,3,1);
%     % b.plot(b.Z .* b.issurface + 0.10); colorbar off; title("Surfaces"+newline+"(not edge, ground, or tiny)");
%     % ax1.FontSize = fontsize; 
%     % hold on;
%     % surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
%     % hold off;
%     % view(3); camlight('headlight'); view(2);
%     ax2 = subplot(1,2,1);
%     % b.plot(b.Z .* b.isroof + 0.10); colorbar off; title("Roofs"+newline+"(surfaces that can fit >1 SPV)");
%     b.plot(b.Z .* b.isroof + 0.10); colorbar off; title("Roofs"+newline+"(not edge, ground, or tiny)");
%     hold on;
%     surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
%     ax2.FontSize = fontsize;
%     hold off;
%     view(3); camlight('headlight'); view(2);
%     ax3 = subplot(1,2,2);
%     b.plot(b.Z .* b.isroof_viable + 0.10); colorbar off; title("Viable roofs"+newline+"(roofs that are flat/slanted)");
%     hold on;
%     surf(b.X,b.Y,b.Z,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
%     ax3.FontSize = fontsize;
%     hold off;
%     view(3); camlight('headlight'); view(2);
%     width=1920;height=1080;
%     set(gcf,'position',[0,0,width,height]);
%     saveas(figure(1),folder+"comparisons/"+name+"_identified.png");
%     close(1);
    
    % Roof types, labeled
    figure(1);
    ax1 = subplot(1,2,1); imagesc(b.map_r); axis image;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    title(name+' Satellite Photo');
    ax1.FontSize = fontsize;   
    
    
    rooftype = zeros(size(b.roofs)); % Will be labeled according to type
    for i = 1:max(b.roofs,[],'all')
        type = b.info.type(b.info.roof==i);
        rooftype(b.roofs==i) = 3*(type=='flat') + 2*(type=='slanted') + 1*(type=='irregular');
    end
    ax2 = subplot(1,2,2); imagesc(rooftype); axis image;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    % cmap = cat(1,0.15*[1,1,1],pink(3));
    cmap = cat(1,0.15*[1,1,1],[0.5774,0,0],[0.8165    0.7165    0.3774],[0.9,0.9,0.8]);
    colormap(cmap);
    L = line(ones(size(cmap,1)-1),ones(size(cmap,1)-1),'LineWidth',10);
    set(L,{'color'},mat2cell(cmap(2:end,:),ones(1,size(cmap,1)-1),size(cmap,2)));            % set the colors according to cmap
    legend(L, {'irregular', 'slanted','flat'}, 'Location','southwest')
    title(name+" Rooftop Types");
    ax = gca; ax.FontSize = fontsize;
    
    width=1920;height=1080;
    set(gcf,'position',[0,0,width,height]);
    saveas(figure(1),folder+"maps/"+name+"_roof_types.png");
    close(1);
end