function [Z_terrain] = ground_gen(Z,X,Y,pix2m)
%GROUND_GEN Creates a 2D approximation of the terrain, by detecting patches
%   of "ground" between buildings, interpolating to connect them, and
%   fitting a 2D surface to the result.
    
    % Goal: Apply edge detection to height map, then select large, low-lying 
    % blobs to be the "ground".
    m_max = 5/12;
    panel_size_pix = 40;
    Z_faketerrain = ones(size(Z)) * -99999;
    issurface = find_surfaces(Z,Z_faketerrain,pix2m, m_max, panel_size_pix);
    issurface = bwareaopen(issurface, ceil(numel(Z)*0.0025),4); % Ignore every blob that takes <0.25% of entire image
    surfaces = bwlabel(issurface,4);
    n_surfaces = max(surfaces,[],'all');

    % Loop through blobs and mark low-lying ones.
    slice_height_low = ground_maxheight() + 1.0;
    i_low = zeros(1,n_surfaces)+9999;   % Preallocating. Hopefully there won't actually be 9999 slices.
    for i = 1:n_surfaces
        slice = double(surfaces==i);
        slice(slice==0) = nan;
    %     slice_height = nanmax(Z .* slice,[],'all');
    %     slice_height = nanmean(Z .* slice,'all');
        slice_height = prctile(Z .* slice,90,'all');
        if slice_height < slice_height_low
            i_low(i) = i;
%             disp(slice_height);
        end
    end

    % Fit a plane to the surface at (surfaces==i_lowest)!
    isground = double(ismember(surfaces,i_low));
    isground(isground==0) = nan;
    Z_ground = Z .* isground;

    % Resize for performance
    smol_scale = 0.2;
    smol_mode = 'nearest';
    X_smol = imresize(X,smol_scale,smol_mode);
    Y_smol = imresize(Y,smol_scale,smol_mode);
    Z_ground_smol = imresize(Z_ground,smol_scale,smol_mode);

    % Convert X, Y, Z to 1D and remove NaNs
    x = X(:);   y = Y(:);   z = Z_ground(:);
    x = x(~isnan(z)); y = y(~isnan(z)); z = z(~isnan(z));
    % Convert X, Y, Z to 1D and remove NaNs
    x_smol = X_smol(:);   y_smol = Y_smol(:);   z_smol = Z_ground_smol(:);
    x_smol = x_smol(~isnan(z_smol)); y_smol = y_smol(~isnan(z_smol)); z_smol = z_smol(~isnan(z_smol));
    
    % Fit!
    % sf = fit([x, y], z, 'linearinterp');
    % sf = fit([x, y], sf(x,y), 'poly22');
    sf_smol = fit([x_smol,y_smol],z_smol,'linearinterp');
    sf_smol = fit([x_smol,y_smol],sf_smol(x_smol,y_smol),'poly22');
    Z_terrain = sf_smol(X,Y) + (Z*0);
    Z_terrain_lims = [min(Z_terrain,[],'all'), max(Z_terrain,[],'all')];
    
    % (Plots available in the `ground_gen_...._test.m` functions. Implement
    % those when needed!)
%     figure(1);
%     fontsize = 15;
%     ax1 = subplot(1,2,1); imagesc(Z_ground_smol); axis image; colorbar; caxis(Z_terrain_lims);
%     title('Detected Ground Height (m)');
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
% %     ax2 = subplot(1,2,2); surf(X,Y,Z, 'EdgeAlpha','0.0','FaceAlpha','0.5'); daspect([1 1 1]); set(gca,'ydir','reverse');
% %     hold on; surf(X,Y,Z_terrain); hold off;
% %     title('Building with Terrain Height Map (m)');
%     ax2 = subplot(1,2,2); imagesc(Z_terrain); axis image; colorbar; caxis(Z_terrain_lims);
%     title('Projected Terrain Height (m)');
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     ax1.FontSize = fontsize;
%     ax2.FontSize = fontsize;

    function maxheight = ground_maxheight()
        % Goal: Apply "min-within-radius" filter (i.e. erosion) to the 
        %   Z-map, and use this to find the maximum height of the ground 
        %   (hopefully).

        % Create a strel of 10m radius, and erode!
        se = strel('disk',ceil(10/pix2m));
        Z_eroded = Z;
        Z_eroded(isnan(Z_eroded)) = 99999;  % Remove NaNs before eroding.
        Z_eroded = imerode(Z_eroded,se);
        Z_eroded = Z_eroded + 0*Z;          % Put NaNs back in.
        
        % Edge detection?
        m_max = 12/12;
        isblob = find_surfaces(Z_eroded,Z_faketerrain,pix2m, m_max, panel_size_pix);
        isblob = double(bwareaopen(isblob, ceil(numel(Z)*0.05),4)); % Ignore every blob that takes <5% of entire image
        Z_blobs = Z_eroded .* isblob;
        Z_blobs(Z_blobs == 0) = nan;      % Edges and deleted blobs are masked out

        % Loop through blobs and directly measure heights around them! If lower
        % than blob height, then delete the blob since it's probably not the
        % ground.
        blobs = bwlabel(isblob,4);
        n_blobs = max(blobs,[],'all');
        for i = 1:n_blobs
            slice = (blobs==i);
            slice_surround = double(imdilate(slice,strel('disk',10)) - slice);  % Region 10 pixels around slice
            % Convert `slice` and `slice_surround` to proper masks
            slice_mask = double(slice);
            slice_mask(slice_mask==0) = nan;
            slice_surround(slice_surround==0) = nan;
            height_surround = nanmean(slice_surround .* Z_blobs,'all');
            % If slice is significantly higher than surroundings, then
            % remove it from consideration. (If height_surround is nan (i.e. if
            % there's only one blob), then this condition returns False and assumes
            % the blob is the ground.)
            if nanmean(slice_mask .* Z_blobs,'all') > (height_surround + 1) 
                blobs = blobs - slice*i;    
            end
        end
        isblob = double(logical(blobs));
        isblob(isblob==0) = nan;
        maxheight = prctile(Z_blobs .* isblob, 99, 'all');
        
%         % Plot
%         figure(2);
%         subplot(2,3,1); imagesc(Z); axis image; colorbar;
%         subplot(2,3,2); imagesc(Z_eroded, [nanmin(Z,[],'all'), nanmax(Z,[],'all')]); axis image; colorbar;
%         subplot(2,3,3); imagesc(isblob); axis image; colorbar;
%         subplot(2,3,4); imagesc(Z_blobs .* isblob); axis image; colorbar;
% %         subplot(2,3,5:6); plot(x,v,x,v2,'r--',x,GROUNDMAX+x*0); title("Percentiles");
    end

end

