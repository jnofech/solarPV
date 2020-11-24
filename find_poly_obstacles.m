function [hole_corners] = find_poly_obstacles(roof,corners,pix2m)
%FIND_POLY_OBSTACLES() approximates polygon shapes for 
%   obstacles on a rooftop.
% 
%   Parameters:
%   -----------
%   roof : double (2D array)
%       Binary image indicating roof.
%   corners : double (2D array)
%       An n-by-2 array (where n is the number of
%       "vertices" on the polygon), in order of
%       connectivity.
%       This is for identifying obstacles that are not
%       properly enclosed within the rooftop.
%
%   Returns:
%   --------
%   hole_corners : cell (of 2D arrays)
%       Contains n-by-2 arrays (where n is the number of
%       "vertices" on each polygon), in order of
%       connectivity. Each array corresponds to an
%       individual hole.

            % Create solid binary out of polygon. Erode by about 0.5m(?).
            % roof_ishole = ~roof .* roof_poly_eroded;
            % Remove tiny holes
            % roof_holes = bwlabel(roof_ishole,4);

            % Loop through holes; find polygons for each hole.


    % Draw rooftop polygon!
    map_polygon = zeros(size(roof));
    for j = 1:size(corners,1)
        x_current = corners(j,1);
        y_current = corners(j,2);
        if j~=size(corners,1)
            x_next = corners(j+1,1);
            y_next = corners(j+1,2);
        else
            x_next = corners(1,1);
            y_next = corners(1,2);
        end
        % Draw a line from current corner to next corner!
        [x_b,y_b] = bresenham(x_current,y_current,x_next,y_next);
        line_len = length(x_b);
        for k = 1:line_len
            x = x_b(k);
            y = y_b(k);
            map_polygon(y,x) = 1;
        end
    end

    % Create solid binary out of polygon. Erode by about 0.5m.
    roof_poly = imfill(map_polygon,'holes');
    roof_poly = imerode(roof_poly,strel('disk',ceil(0.5/pix2m)));
    roof_ishole = roof_poly .* ~roof;

    % Identify individual obstacles ("holes" in image)!
    roof_holes = bwlabel(roof_ishole,4);
    n_holes = max(roof_holes,[],'all');
    
    % Loop through obstacles; find polygon corner locations.
    hole_corners  = cell(n_holes,1);   % Will hold the 1x2 outputs `corners` of each hole.
    for ii = 1:n_holes
        corners = pgonCorners(roof_holes==ii,10);
        corners = circshift(corners,1,2);   % Swaps `x` and `y` columns
        hole_corners{ii} = corners;
    end
    
    % Force obstacle polygons to be clockwise, if this was not already the case
    for ii = 1:size(hole_corners,1)
        % Check each obstacle (`hole_corners{i}`) individually.
        if ~is_clockwise(hole_corners{ii})
            hole_corners{ii} = flipud(hole_corners{ii});
        end
    end
end