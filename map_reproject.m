function varargout = map_reproject(map,map_base,X1,X2,Y1,Y2,X1b,X2b,Y1b,Y2b)
%MAP_REPROJECT reprojects a 2D or 3D array `map` to the orientation of the
%   2D or 3D array `map_base`, given two pairs of corresponding
%   coordinates.
%   Assumes that the image in `map` is rotated from `map_base` by some
%   multiple of 90 degrees. Other angles (e.g. an image of `map_base`
%   canted by 45 degrees) will not work.
%
%   Parameters:
%   -----------
%   map : double (2D or 3D array)
%       Map to be reprojected.
%   map_base : double (2D or 3D mapay)
%       Map to which `map` is reprojected.
%   X1,X2,Y1,Y2 : double
%       Pixel coordinates of two specific points on `map`.
%   X1b,X2b,Y1b,Y2b : double
%       Pixel coordinates of these same points on `map_base`.
%
%   Outputs:
%   --------
%   [map_r] : 
%       Reprojected version of `map`. There will likely be a lot of
%       `nan` values.
%   [map_r, X1_r, X2_r, Y1_r, Y2_r] : 
%       Includes coordinates for the two points on `map_r`.

    [ymax_b, xmax_b, zmax_b] = size(map_base);
    
    % Rotate `map`
    coords_direction   = string(""+double(X1<X2)+double(Y1<Y2));  % 1 indicates "increasing".
    coords_direction_b = string(""+double(X1b<X2b)+double(Y1b<Y2b));
    coords_sequence = ["00"; "10"; "11"; "01"];             % Sequence of `coords_direction` in steps of 90deg CW
    rotations_cw = find(coords_sequence == coords_direction_b) - find(coords_sequence == coords_direction);
    rotation_CCW = rotations_cw * -90;
    % ^ Number of CW rotations to reproject cmap to Zmap
    [map,coord1,coord2] = map_rotate(map,rotation_CCW,[X1,Y1],[X2,Y2]);
    X1 = coord1(1); Y1 = coord1(2);
    X2 = coord2(1); Y2 = coord2(2);
    
    
    [ymax, xmax, zmax] = size(map);
    
    % X-direction: Set new bounds, and shift `map`,`X1`,`X2` across these new bounds.
    X1_old = X1; X2_old = X2; xmin_old = 1; xmax_old = xmax;    % BACKUP info for `map`.
    xmin = 1;                                                   % Min x-value of `map_r`.
    xmax = abs(X2_old-X1_old)/abs(X2b-X1b)*(xmax_b-1) + 1;      % Max x-value of `map_r`.
    X1 = (X1b-1)/(xmax_b-1)*(xmax-1)+1;         % x-value of first coordinate in `map_r`.
    X2 = X1 + (X2-X1_old);                      % x-value of second coordinate in `map_r`.
    
    % Slide `map` into `map_r`, and prepare to crop where necessary.
    X_diff = X1 - X1_old;           % The X-shift from old colour map to new colour map
    x_range = xmin_old:xmax_old;    % Range of x-values in `map`
    x_allowed = find(x_range+X_diff>=xmin & x_range+X_diff<=xmax);  % x-values in `map` that fit within `map_r`'s bounds
    x_lowbound = x_allowed(1);      % The bounds, in `map`
    x_highbound = x_allowed(end);
%     disp("BEFORE: (xmin, xmax) = ("+xmin_old+", "+xmax_old+"), and (X_1, X_2) = "+X_1_old+", "+X_2_old+").");
    % NOTE: Difference between x_highbound & x_lowbound MAY NOT EQUAL difference
    % between xmax and xmin. This is because X_diff is not an integer, and may
    % "shift" x_range to its `map_r` version such that an extra integer is 
    % squeezed in between (or unable to fit in between) xmax and xmin
    % (which are also not integers).
    
    % Y-direction: Exact same thing as above.
    Y1_old = Y1; Y2_old = Y2; ymin_old = 1; ymax_old = ymax;
    ymin = 1;
    ymax = abs(Y2-Y1)/abs(Y2b-Y1b)*(ymax_b-1) + 1;
    Y1 = (Y1b-1)/(ymax_b-1)*(ymax-1)+1;
    Y2 = Y1 + (Y2-Y1_old);
    %
    Y_diff = Y1 - Y1_old;
    y_range = ymin_old:ymax_old;
    y_allowed = find(y_range+Y_diff>=ymin & y_range+Y_diff<=ymax);
    y_lowbound = y_allowed(1);
    y_highbound = y_allowed(end);
%     disp("BEFORE: (ymin, ymax) = ("+ymin_old+", "+ymax_old+"), and (Y1, Y2) = "+Y1_old+", "+Y2_old+").");

    % Crop the array accordingly!
%     disp("x low/high bounds in `map`: "+x_lowbound+", "+x_highbound);
%     disp("y low/high bounds in `map`: "+y_lowbound+", "+y_highbound);
    map = imcrop(map,[x_lowbound, y_lowbound, x_highbound-x_lowbound, y_highbound-y_lowbound]);
%     disp("Size of `map` after cropping, but before reprojection: "+size(map,1)+", "+size(map,2)+", "+size(map,3)+".");


    % Create `map_r`.
    xmax = round(xmax);             % Do the rounding AFTER new coords have been accurately calculated.
    ymax = round(ymax);             % Do the rounding AFTER new coords have been accurately calculated.
    map_r = nan(ymax,xmax,zmax);
    map_r(round(1+Y_diff+y_lowbound-ymin_old) : (round(1+Y_diff+y_lowbound-ymin_old)+size(map,1)-1), ...
          round(1+X_diff+x_lowbound-xmin_old) : (round(1+X_diff+x_lowbound-xmin_old)+size(map,2)-1), ...
           : ) = map;               % Fill `map_r` with `map`'s values!
    map_r = uint8(map_r);
%     disp("AFTER: (xmin, xmax) = ("+xmin+", "+xmax+"), and (X1, X2) = "+X1+", "+X2+").");
%     disp("AFTER: (ymin, ymax) = ("+ymin+", "+ymax+"), and (Y1, Y2) = "+Y1+", "+Y2+").");
       
       
    % Return results  
    switch nargout        
        case 1
            varargout{1} = map_r;
        case 2
            varargout{1} = map_r;
            varargout{2} = rotation_CCW;
        case 5
            varargout{1} = map_r;
            varargout{2} = X1;
            varargout{3} = X2;
            varargout{4} = Y1;
            varargout{5} = Y2;
        case 6
            varargout{1} = map_r;
            varargout{2} = X1;
            varargout{3} = X2;
            varargout{4} = Y1;
            varargout{5} = Y2;
            varargout{6} = rotation_CCW;
        otherwise
            error("Can output `[map_r]` or `[map_r, X1_r, X2_r, Y1_r, Y2_r]`.");
    end
        
end

