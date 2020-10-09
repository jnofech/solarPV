function [issurface,isroof] = find_surfaces(Z,Z_terrain,pix2m,m_max,minsize_pix,roof_mindiameter)
%FIND_SURFACES() identifies surfaces from a height map; i.e. any structure 
%   that is not the ground, an edge, or tiny.
%
%   Parameters:
%   -----------
%   Z : double (2D array)
%       Height map of the building.
%   Z_terrain : double (2D array)
%       Height map of the terrain.
%   pix2m : double
%       Conversion factor between pixels and metres for x-, y-, and z-axes.
%   m_max : double
%       Slope identifying "steep" regions. Narrow "steep" regions are
%       treated as edges.
%       Recommended: 12/12.
%   minsize_pix : double
%       Minimum number of pixels that a usable "surface" must have (e.g.
%       minimum area of a solar panel or rooftop).
%
%   Returns:
%   --------
%   issurface : logical (2D array)
%       Indicates any structure that is not the ground, an edge, or tiny.

    % Initialization
    u = symunit;
    % Find gradients from height map (NOT actual face slopes)
    [Zx,Zy] = gradient(Z, pix2m);     % Finds gradient, with each pixel representing `pix2m_h` metres.
    m = sqrt( Zx.^2 + Zy.^2 );        % Magnitude of the slope (i.e. slope in direction of deepest ascent)
    
    % Edge detection
    Z_steep = (m > m_max);  % 1 where slope is too STEEP.
        % Quickly clean!
        Z_steep = bwareaopen(Z_steep,ceil(numel(Z)*1.9e-5));   % Quickly cleaning "noise", where slope is unintentionally high.
    se_1pixel = strel('disk',1);
    Z_steepsurface = imerode(Z_steep,se_1pixel);   % Steep surfaces, NOT edges.
    Z_edge = Z_steep - Z_steepsurface;
    Z_edge = imdilate(Z_edge,se_1pixel);           % Increases edge thickness by 1 on either side.

    % Create `Z_surface`-- Identify every non-ground surface
    Z_isground = Z < (Z_terrain + 3.0);
    issurface = double(~Z_isground & ~Z_edge & ~isnan(Z));     % Any non-edge surfaces, regardless of steepness and lumpiness.

    % Remove structures with less area than minimum roof area.
    issurface = bwareaopen(issurface,minsize_pix);   % Any surface that exceeds specified area.
    
    % Roof detection! (Surfaces that are large enough to fit a tiny solar
    % panel under California regulations regarding edge distance)
    roofs_temp = blob_clean(issurface,roof_mindiameter);
    isroof = logical(roofs_temp);
    
%     % Plot me a thing
%     figure(1);
%     fontsize = 15;
%     ax1 = subplot(1,2,1); imagesc(Z_edge); axis image; colorbar;
%     title('Edges');
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     ax1.FontSize = fontsize;
%     ax2 = subplot(1,2,2); imagesc(Z_edge*0); axis image; colorbar;  % Replace with b.map
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     ax2.FontSize = fontsize;    
end

