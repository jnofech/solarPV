function [shaded] = shading_raytrace(b,alt,azi,scale)
% Given a Building object (b) and the altitude/azimuth towards the Sun,
%   returns a binary map indicating whether a pixel is illuminated (false) 
%   or shaded (true). Can scale the heightmap for performance.

    % (Choose ONE!)
        % Heightmap from Building object
        Z = b.Z_nomask;                     % Height map, in metres.
        Z = Z / b.pix2m;                    % Height map, in pixel units.
%         % Heightmap, inputted directly
%         Z = b;

    % Heightmap resizing
    Z_backup = Z;
    if scale < 1
        Z = imresize(Z,scale,'nearest');
    elseif scale == 1
        % do nothing
    else
        % also do nothing
        warning('shading_raytrace.m  : `scale` should be <=1.');
    end
    
    xaxis = 1:size(Z,2);
    yaxis = 1:size(Z,1);
    shaded = zeros(size(Z),'logical');  % (Will be) output binary array.

    % Get ray vector.
    % (dx,dy,dz represent the Cartesian unit vector pointing towards the
    % sun.)
    comp = 1*exp(deg2rad(-90-azi)*1i);  % Complex number, where 0deg=E, 90deg=N, etc
    dx = real(comp);
    dy = -imag(comp);                   % Flip to account for MATLAB's method
    dz = 1*sind(alt);

    % Create bresenham line.
    % (Essentially a pixelated line connecting two points in an image.)
    x0=0;   y0=0;
    maxdist = norm(size(Z));  % Maximum pixel distance of two points in the
                              % height map.
    x1 = (x0+maxdist*dx);   % Arbitrarily far-away x-coordinate in Sun's direction.
    y1 = (y0+maxdist*dy);   % Arbitrarily far-away y-coordinate in Sun's direction.
    [x_b,y_b] = bresenham(x0,y0,x1,y1);     % Bresenham line in Sun's direction.
    
    % Get dz/dx, dz/dy of ray direction
    zx = dz/dx;     zy = dz/dy;
    if ~isfinite(zx)
        zx = 0;     % Prevents error if ray has no East/West component
    end
    if ~isfinite(zy)
        zy = 0;     % Prevents error if ray has no North/South component
    end
    
    % Convert bresenham-related vars to int16, for quicker calculations.
    xaxis   = int16(xaxis);     yaxis   = int16(yaxis);
    x_b = int16(x_b);   y_b = int16(y_b);
    x_min = int16(0);  x_max = int16(size(Z,2));
    y_min = int16(0);  y_max = int16(size(Z,1));
    
    % Initialize index cells
    idx_xcell = cell(x_max,1);
    idx_ycell = cell(y_max,1);
    for x0 = xaxis
        x_bp = x_b + x0;   	% Get Bresenham line's x-values for each x-position in heightmap.
        idx_xcell{x0} = (x_bp > x_min) & (x_bp <= x_max);   % Track where Bres. line falls outside heightmap's x-range.
    end
    for y0 = yaxis
        y_bp = y_b + y0;   	% Get Bresenham line's y-values for each y-position in heightmap.
        idx_ycell{y0} = (y_bp > y_min) & (y_bp <= y_max);   % Track where Bres. line falls outside heightmap's y-range.
    end
    
    % Loop through all pixels in heightmap!
    if ~(abs(dx) == 0 && abs(dy) ==0)
        for x0 = xaxis
            x_bp = x_b + x0;   	% Bresenham x-values at each x-position.
            for y0 = yaxis
                y_bp = y_b + y0;    % Bresenham y-values at each y-position.

                % Combine x- and y- "good" indices
                % (i.e. indices where the Bresenham line is contained
                % within heightmap).
                index = idx_xcell{x0} & idx_ycell{y0};
                x_bp = x_bp(index);
                y_bp = y_bp(index);

                % Get "ray height" z_b along bressenham line; z = z0+ (dz/dx)(x-x0) + (dz/dy)(y-y1)
                z0 = Z(y0,x0);
                % Projected ray heights (relative)
                z_b = z0 + zx*double(x_bp-x0) + zy*double(y_bp-y0);
                % Get ACTUAL heights along line
                z_a = Z(sub2ind(size(Z),y_bp,x_bp));      % Building heights

                % Compare the two heights-- does the ray intersect with anything?
                if any(z_a > z_b)
                    shaded(y0,x0) = true;
                else
                    shaded(y0,x0) = false;
                end
            end
        end
    elseif abs(dx) == 0 && abs(dy) ==0 && dz>0
        % entire map is illuminated
        % OLD: 1
        shaded(:,:) = false;
    elseif abs(dx) == 0 && abs(dy) ==0 && dz<=0
        % entire map is dark 
        % OLD: 0
        shaded(:,:) = true;
    end
    
    % Resize "shaded" to original dimensions!
    shaded = imresize(shaded, [size(Z_backup,1),size(Z_backup,2)]);
end