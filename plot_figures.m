function plot_figures(b,pix_length,lat,output_path,output_path_save)   
    % Initialize
    fontsize = 14; 
    name = b.name;
    % Establish paths
    output_path = output_path+"figures\";
    temp_paths = [output_path+"maps\";
                  output_path+"identified\";
                  output_path+"polygons\"+name+"\";
                  output_path+"comparisons\"];
    for i = 1:size(temp_paths,1)
        temp_path = temp_paths(i);
        if ~exist(temp_path,'dir')
            mkdir(temp_path);
        end
    end
    
    % ~~~ Begin plotting! ~~~

    % Height map, with mask
    figure(1);
    b.plot(b.Z);
    ax = gca; ax.FontSize = fontsize;
    saveas(figure(1),output_path+"maps\"+name+"_hmap_masked.png");

    % Height map
    figure(1);
    b.plot(b.Z_nomask);
    ax = gca; ax.FontSize = fontsize;
    saveas(figure(1),output_path+"maps\"+name+"_hmap.png");

    % Height+Terrain map
    hold on; alpha 0.5; surf(b.X,b.Y,b.terrain_nomask); hold off;
    ax = gca; ax.FontSize = fontsize;
    view(3);
    title("Height Map of "+name+" (m),"+newline+"with approximate terrain (black)");
    saveas(figure(1),output_path+"maps/"+name+"_hmap_terrain.png");
    close(1);

    % Colour map
    figure(1);
    b.plot(b.map_r);
    ax = gca; ax.FontSize = fontsize;
    saveas(figure(1),output_path+"maps/"+name+"_cmap.png");
    close(1);

    % Height+Colour maps combined
    figure(1);
    width=1920;height=1080;
    set(gcf,'position',[0,0,width,height]);
    subplot(1,2,1);
    b.plot(b.Z_nomask);
    ax = gca; ax.FontSize = fontsize;
    subplot(1,2,2);
    b.plot(b.map_r);
    ax = gca; ax.FontSize = fontsize;
    saveas(figure(1),output_path+"maps/"+name+"_both.png");
    close(1);
    
    % Height+Colour maps combined, masked
    figure(1);
    width=1920;height=1080;
    set(gcf,'position',[0,0,width,height]);
    subplot(1,2,1);
    b.plot(b.Z);
    ax = gca; ax.FontSize = fontsize;
    subplot(1,2,2);
    b.plot(b.map_r);
    ax = gca; ax.FontSize = fontsize;
    saveas(figure(1),output_path+"maps/"+name+"_both_masked.png");
    close(1);

    % % Surfaces
    % figure(1);
    % b.plot(b.Z .* b.issurface);
    % ax = gca; ax.FontSize = fontsize;
    % title("Surfaces of "+name);
    % saveas(figure(1),folder+"identified/"+name+"_surfaces.png");
    % close(1);

    % Roofs
    figure(1);
    b.plot(b.Z .* b.isroof_viable);
    ax = gca; ax.FontSize = fontsize;
    title(name+": Panel-Viable Rooftops");
    saveas(figure(1),output_path+"identified/"+name+"_roofs.png");
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
    % b.plot(b.Z_nomask .* b.isroof + 0.10); colorbar off; title("Roofs"+newline+"(surfaces that can fit >1 SPV)");
    b.plot(b.Z_nomask .* b.isroof + 0.10); colorbar off; title(name+": All Roofs"+newline+"(anything that isn't an edge, the ground, or tiny)");
    hold on;
    surf(b.X,b.Y,b.Z_nomask,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
    ax2.FontSize = fontsize;
    hold off;
    view(3); camlight('headlight'); view(2);
    ax3 = subplot(1,2,2);
    b.plot(b.Z_nomask .* b.isroof_viable + 0.10); colorbar off; title(name+": Viable roofs"+newline+"(roofs that are flat/slanted)");
    hold on;
    surf(b.X,b.Y,b.Z_nomask,'edgecolor','none','facecolor',[0.4 0.4 0.4]); set(gca,'ydir','reverse'); daspect([1 1 1]); 
    ax3.FontSize = fontsize;
    hold off;
    view(3); camlight('headlight'); view(2);
    width=1920;height=1080;
    set(gcf,'position',[0,0,width,height]);
    saveas(figure(1),output_path+"comparisons/"+name+"_identified.png");
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
    colormap(cmap); caxis([0,3]);
    L = line(ones(size(cmap,1)-1),ones(size(cmap,1)-1),'LineWidth',10);
    set(L,{'color'},mat2cell(cmap(2:end,:),ones(1,size(cmap,1)-1),size(cmap,2)));            % set the colors according to cmap
    legend(L, {'irregular', 'slanted','flat'}, 'Location','southwest')
    ax = gca; ax.FontSize = fontsize*2/3;
    title(name+": Rooftop Classifications"+newline+"('irregular' roofs are discarded)");
    saveas(figure(1),output_path+"maps/"+name+"_roof_types.png");
    close(1);
    
    % Roof types, labeled, and compared with satellite photos
    
            % Height+Colour maps combined
            figure(1);
            subplot(1,2,1);
            b.plot(b.Z);
            ax = gca; ax.FontSize = fontsize;
            subplot(1,2,2);
%             b.plot(b.map_r);
%             ax = gca; ax.FontSize = fontsize;
%             saveas(figure(1),output_path+"maps/"+name+"_both.png");
%             close(1);
    
    figure(1);
    width=1920;height=1080;
    set(gcf,'position',[0,0,width,height]);
    subplot(1,2,1);
    b.plot(b.map_r);
    ax = gca; ax.FontSize = fontsize;
    subplot(1,2,2);    
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
    title(name+": Rooftop Classifications"+newline+"('irregular' roofs are discarded)");
    saveas(figure(1),output_path+"comparisons/"+name+"_cmap_vs_roof_types.png");
    close(1);


    % Shadow study, full year
%     lat = 53.5461;        % Edmonton's latitude (deg, North)
    % lon = 113.4938;       % Edmonton's longitude (deg, West)
    day = NaN;            % 172 for summer solstice (Jun 21), 355 for winter solstice (Dec 21), 81 for equinox (Mar 21, Sep 21)
    dates = [81, 172, 264, 355];     % Mar 21, Jun 21, Sep 21, Dec 21 (Mar 21, Sep 21 are identical)
    hours = 1:24;               % Hourly checks
    load(output_path_save+name+"_shadows_"+pix_length+".mat");
    % Get hillshade
    for i = 1:numel(dates)
        day = dates(i);
        dec = 23.45*sind(360/365*(284+day));    % (approximate) Solar declination; ranging from 23.45 to -23.45
        H_angle_rise = -acosd( -tand(lat) * tand(dec) );    % Hour angle at sunrise
        hour_rise = H_angle_rise / 60 / 0.25 + 12.00;       % Hour at sunrise
        for j = 1:numel(hours)
            hour_AST = hours(j);
            if (day == 172) && (hour_AST == 12)
                % Get solar angles
                H_angle = 0.25*(hour_AST - 12.00)*60;   % Hour angle; i.e. 0.25deg per minute since Noon
                dec = 23.45*sind(360/365*(284+day+hour_AST/24));    % Solar declination; ranging from 23.45 to -23.45
                alt_noon = 90-lat+dec;                  % Altitude of sun in sky AT SOLAR NOON
                alt = asind( cosd(lat)*cosd(dec).*cosd(H_angle) + sind(lat)*sind(dec) ); % Altitude
                azi = asind( sind(H_angle).*cosd(dec)/cosd(alt) );                       % Azimuth (deg CW from South)
                % Get hillshade
                hill = hillshade(b.Z_nomask,b.X(1,:),b.Y(:,1),'azimuth',azi+180,'altitude',alt)/255;    % Input azi is in deg CW from North.
            end
        end
    end
    % Plot shadows with hillshade (June 21 at noon is used for hillshade)
    figure(1);
    imshow(hill);
    color = [41, 93, 204] / 255;
    color_im = zeros([size(b.Z),3]);
    color_im(:,:,1) = color(1);
    color_im(:,:,2) = color(2);
    color_im(:,:,3) = color(3);
    shaded_final = sum(shadows,3);
    shaded_final = shaded_final / size(shadows,3);      % 0 indicates light, 1 indicates darkness
    shaded_alpha = (shaded_final - min(min(shaded_final,[],'all'),0.5)) / (1 - 0.5);     % Display transparency for shadows.
    hold on; im_shaded = imshow(color_im); hold off;
    set(im_shaded, 'AlphaData', shaded_alpha);
    title(name+": Full-year shadow study");
    saveas(figure(1),output_path+"maps/"+name+"_shadows_fullyear.png");
    close(1);

    % Shadow study, one day at a time
%     lat = 53.5461;        % Edmonton's latitude (deg, North)
    % lon = 113.4938;       % Edmonton's longitude (deg, West)
    day = NaN;            % 172 for summer solstice (Jun 21), 355 for winter solstice (Dec 21), 81 for equinox (Mar 21, Sep 21)
    dates = [81, 172, 264, 355];     % Mar 21, Jun 21, Sep 21, Dec 21 (Mar 21, Sep 21 are identical)
    hours = 1:24;               % Hourly checks
    % One date at a time
    for date = dates
        % Shadow study!
        load(output_path_save+name+"_shadows_"+pix_length+".mat");
        shaded_final = zeros(size(b.Z));

        for i = 1:numel(dates)
            day = dates(i);
            dec = 23.45*sind(360/365*(284+day));    % (approximate) Solar declination; ranging from 23.45 to -23.45
            H_angle_rise = -acosd( -tand(lat) * tand(dec) );    % Hour angle at sunrise
            hour_rise = H_angle_rise / 60 / 0.25 + 12.00;       % Hour at sunrise

            for j = 1:numel(hours)
                hour_AST = hours(j);

                if (day == date)
                    % Get solar angles
                    H_angle = 0.25*(hour_AST - 12.00)*60;   % Hour angle; i.e. 0.25deg per minute since Noon
                    dec = 23.45*sind(360/365*(284+day+hour_AST/24));    % Solar declination; ranging from 23.45 to -23.45
                    alt_noon = 90-lat+dec;                  % Altitude of sun in sky AT SOLAR NOON
                    alt = asind( cosd(lat)*cosd(dec).*cosd(H_angle) + sind(lat)*sind(dec) ); % Altitude
                    azi = asind( sind(H_angle).*cosd(dec)/cosd(alt) );                       % Azimuth (deg CW from South)

                    % Get shadow
                    shaded = double(shadows(:,:,(i-1)*numel(hours)+j));
                    shaded_final = shaded_final + shaded;

                    if hour_AST == 12
                        % Get hillshade at noon
                        hill = hillshade(b.Z_nomask,b.X(1,:),b.Y(:,1),'azimuth',azi+180,'altitude',alt)/255;    % Input azi is in deg CW from North.
                    end
                end
            end
        end

        % Plot shadows with (noon) hillshades
        figure(1);
        imshow(hill);
        color = [41, 93, 204] / 255;
        color_im = zeros([size(shaded),3]);
        color_im(:,:,1) = color(1);
        color_im(:,:,2) = color(2);
        color_im(:,:,3) = color(3);
        shaded_final = shaded_final / max(shaded_final,[],'all');
        shaded_alpha = (shaded_final - min(min(shaded_final,[],'all'),0.5)) / (max(shaded_final,[],'all') - min(min(shaded_final,[],'all'),0.5));
        hold on; im_shaded = imshow(color_im); hold off;
        set(im_shaded, 'AlphaData', shaded_alpha);
        title(name+": Shadows throughout day "+date+" of 365");
        saveas(figure(1),output_path+"maps/"+name+"_shadows_"+date+".png");
        close(1);
    end
    
    % Plot all polygons
    n_roofs = max(b.roofs,[],'all');
    for i = 1:n_roofs
        figure(1);
        width=1366;height=768;
        set(gcf,'position',[200,200,width,height]);

        b.roofplot(i);
        saveas(figure(1),output_path+"polygons\"+name+"\"+name+"_polygon_"+i+".png");
    end
    close(1);
end