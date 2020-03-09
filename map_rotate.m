function varargout = map_rotate(varargin)
%MAP_ROTATE rotates a 2D or 3D array `map`, and any inputted coordinates
%   within these maps, by some specified angle about the Z-axis.
%
%   Parameters:
%   -----------
%   map : double (2D or 3D array)
%       Map to be rotated. Will be facing North afterwards.
%   theta_CCW : double
%       Counterclockwise rotation, in degrees.
%   (OPTIONAL)
%   [x,y] : double (1x2 array)
%       Coordinates (in 2D pixel space) to be rotated as well.
%       Multiple pairs of coordinates can be inputted as separate
%       arrays.
%
%   Outputs:
%   --------
%   [map_rot] : 
%       Rotated version of `map`.
%   [x_rot,y_rot] : 
%       Rotated coordinates. Will return as many coordinate pairs as are
%       inputted.
    if nargin==2
        map = varargin{1};
        theta_CCW = varargin{2};
        x = nan;
        y = nan;
    elseif nargin>=3
        map = varargin{1};
        theta_CCW = varargin{2};
        n_pairs = size(varargin,2)-2;
        x = zeros(1,n_pairs);
        y = zeros(1,n_pairs);
        for i = 1:n_pairs
            x(i) = varargin{2+i}(1);
            y(i) = varargin{2+i}(2);
        end
    else
        error("Please input a 2D or 3D map you would like to see rotated, "...
             +"as well as the number of CCW rotations.");
    end
    
    % 3D rotation matrix
    Rz = [  cosd(theta_CCW), -sind(theta_CCW), 0;
            sind(theta_CCW),  cosd(theta_CCW), 0;
                          0,                0, 1];
    Rz_coord = [  cosd(theta_CCW),  sind(theta_CCW), 0;
                 -sind(theta_CCW),  cosd(theta_CCW), 0;
                                0,                0, 1];    % Y-direction reversed, because of how arrays are set up in MATLAB
    switch size(map,3)
        case 1
            % Assumed to be height map
            Rz       = Rz(1:2,1:2);
            Rz_coord = Rz_coord(1:2,1:2);
            coords = [x;y];
        case 3
            % Assumed to be colour image. Heightmap is used to convert to
            % m, since the colour map isn't perfectly top down.
            coords = [x;y;zeros(size(x))];
        otherwise
            error('For `map` inputs, only (N x M x 1) and (N x M x 3) arrays are supported.')
    end
    
    map_cen = (size(map)/2)';                       % Centre of map
    map_cen([1 2]) = map_cen([2 1]);                % Swap X and Y (since rows, i.e. y's, come first in MATLAB)
    map_rot = imrotate(map,theta_CCW,'bilinear');   % Rotates map, filling gaps with zeros.
    map_rot_cen = (size(map_rot)/2)';               % Centre of rotated map
    map_rot_cen([1 2]) = map_rot_cen([2 1]);                % Swap X and Y (since rows, i.e. y's, come first in MATLAB)
    coords = Rz_coord*(coords-map_cen)+map_rot_cen;

    varargout = cell(1,nargin-1);
    varargout{1} = map_rot;
    if nargin>=3
        for i = 1:n_pairs
            varargout{1+i} = coords(:,i);
        end
    end
end

