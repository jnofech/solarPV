classdef Building
    % BUILDING object, where name is specified as input.
    %   Contains info about 3D mesh and 2D topview of building, as well as
    %   conversions from pixel to real-life units.
    
    properties
        name = 'Null'   % Building name (case-sensitive).
        length;         % Pixel length of 3D map (user-specified).
        mesh;           % 3D mesh info (faces, vertices, normals).
        Z;              % Height map.
        x;              % X-coordinate array, in m.
        y;              % Y-coordinate array, in m.
        X;              % X-coordinate meshgrid.
        Y;              % Y-coordinate meshgrid.
        dZdX;           % Slope in X-direction.
        dZdY;           % Slope in Y-direction.
        m;              % Magnitude of slope.
        pix2m;          % Metres/pixel, for Z-map
        terrain;        % Height map of approximate terrain elevation.
        issurface;      % Binary map indicating non-ground surfaces.
        isroof;         % Binary map indicating roofs.
        isroof_viable   % Binary map indicating roofs that are either flat or slanted (NOT curved or irregular).
        roofs;          % Blob map (from `bwlabel()`) of each roof.
        roof_isbad;     % Regions where each roof (value corresponding to `roofs`) has unreliable values.
        info;           % Table containing information of each blob.
        % Trash pile
        map;            % Colour image of the building.
        map_r;          % Colour image of the building, reprojected to mesh's proportions.
        pix2m_r_x;      % Metres/pixel, for reprojected C-map
        pix2m_r_y;      % Metres/pixel, for reprojected C-map
    end
    
    methods
        function obj = Building(varargin)
            % CONSTRUCT: Allow for variable inputs, as long as they're the
            % building's name and pixel length (i.e. resolution).
            switch nargin
                case 1
                    obj.name = varargin{1};
                    warning("Building(): No pixel length specified; setting to default of 400.");
                    obj.length = 400;
                case 2
                    obj.name = varargin{1};
                    obj.length = varargin{2};
                case 0
                    error("This building has a NAME, you poopnobbler!");
                otherwise
                    error('Use `Building(name,length)`. Length is optional.');
            end
            
            % Now, with the building named: Initialize!
            u = symunit;
            obj.map = map_read(obj.name);
            fname_stl_rot = "output/"+obj.name+"_map_percent2m.mat";
            
            % ~~~~~~~~~~~~~~~ GET HMAP+CMAP ~~~~~~~~~~~~~~~~~~~~
            while true
                if isfile(fname_stl_rot)
                    load(fname_stl_rot,'rotated_CCW');
                    theta_CCW = -rotated_CCW;
                else
                    theta_CCW = 0;
                    rotated_CCW = nan;
                end
                [obj.Z,obj.X,obj.Y,obj.dZdX,obj.dZdY,obj.mesh] = heights_gen(obj.name,obj.length,theta_CCW);
                obj.m = sqrt( (obj.dZdX).^2 + (obj.dZdY).^2 );  % Magnitude of the slope
                obj.x = obj.X(1,:); obj.y = obj.Y(:,1);
                % Get rid of "inf" values
                obj.Z = obj.Z .* isfinite(obj.Z);

                % Convert from pixels to metres! (These arrays will start at
                % ZERO, not 1, since they don't match indices anymore.)
                obj.pix2m = pix_to_m(obj.name,obj.Z,theta_CCW); % Also saves `theta_CCW` into '..._heightmap_..._percent2m.mat'.
                obj.Z = (obj.Z-1)*obj.pix2m;    % Just for consistency, even though heights are all relative to ground anyways.
                obj.X = (obj.X-1)*obj.pix2m;
                obj.Y = (obj.Y-1)*obj.pix2m;
                obj.x = obj.X(1,:); obj.y = obj.Y(:,1);
                [obj.map_r, obj.pix2m_r_x, obj.pix2m_r_y] = pix_to_m(obj.name,obj.map);
                % Pixel conversion factor
                pix = newUnit('pix',obj.pix2m*u.m);
                if theta_CCW == -rotated_CCW
                    break
                end
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Normalize "Z" to height of ground
            obj.terrain = ground_gen(obj.Z,obj.X,obj.Y,obj.pix2m);
%             saveas(figure(1),"output/fancyplot/maps/"+obj.name+"_hmap_terrain.png");
            Z_ground = prctile(obj.terrain,50,'all');
            obj.Z = obj.Z - Z_ground;
            obj.terrain = obj.terrain - Z_ground;
            
            % Solar PV parameters (!!! Maybe set this up in its own
            % function / object? !!!)
            panel_orientation = 'landscape';
            panel_tilt = 0*u.deg;                     % Tilt angle of panel from true horizontal (NOT from rooftop)
            panel_y = 65*u.in;                        % Min. length of rooftop, including regulatory distance from edges
            panel_x = 39*u.in;                        % Min. width of rooftop, including regulatory distance from edges
            switch panel_orientation
                case 'landscape'
                    panel_length = panel_y;
                    panel_width  = panel_x*cos(panel_tilt);
                case 'portrait'
                    panel_length = panel_y*cos(panel_tilt);
                    panel_width  = panel_x;
                otherwise
                    error('`panel_orientation` must be "portrait" or "landscape"');
            end
            % Specify minimum roof dimensions (panel + distance from edges)
            roof_edgedistance = 1.8*u.m;            % Regulatory distance from edge for each PV unit
            panel_length = rewrite(panel_length, u.pix);
            panel_width  = rewrite(panel_width, u.pix);
            roof_edgedistance = rewrite(roof_edgedistance, u.pix);
            roof_minlength = 2*roof_edgedistance+panel_length;
            roof_minwidth  = 2*roof_edgedistance+panel_width; 
%             roof_mindiameter = min([val(roof_minlength),val(roof_minwidth)]);
            roof_mindiameter = min([val(panel_length),val(panel_width)]);
            % Total minimum area of solarPV-feasible rooftop?
            panel_size_pix = ceil(val(panel_length*panel_width));     % Number of pixels in each solar panel.
            roof_minsize_pix = ceil(val(roof_minlength*roof_minwidth));
            % Note: `val(x)` is just `double(separateUnits(x))`.    
            
            % ~~ Roof Identification ~~
            % Identify any region that isn't the ground, an edge, or tiny
            m_max = 12/12;
            obj.issurface = find_surfaces(obj.Z,obj.terrain,obj.pix2m, m_max, roof_minsize_pix);
%             obj.surfaces = bwlabel(obj.issurface,4);
            % Identify any surface that can be treated as a roof
            obj.isroof = find_roofs(obj.issurface, roof_mindiameter);  % Surfaces that may fit an SPV (??? is `blob_clean` doing it correctly?)
            obj.roofs = blob_heights_split(obj.isroof,obj.Z);     % Roofs, properly separated from each other
            obj.isroof = logical(obj.roofs);                      % Binary for the above
            obj.roof_isbad = blob_clean_spikes(obj.isroof, obj.Z); % Roofs, with sudden "spikes" in height cleaned

            % Remove roofs from `surfaces` (re-add after `roof` modifications)
            obj.issurface = obj.issurface - obj.isroof;
            
            % Classify roofs & split some `irregular` roofs!
            [obj.roofs,obj.roof_isbad,obj.info] = blob_classify(obj.roofs,obj.roof_isbad,obj.Z,obj.dZdX,obj.dZdY,obj.pix2m,roof_mindiameter,1);
%             [obj.roofs,obj.roof_isbad,obj.info] = blob_classify_noslopesplit(obj.roofs,obj.roof_isbad,obj.Z,obj.dZdX,obj.dZdY,obj.pix2m,roof_mindiameter,1);
            warning("blob_classify() : stdev threshold for plane fitting arbitrarily set to 0.09, after working with CAM's wavy rooftop. Need to test some more on sloped surfaces and irregular surfaces!");
            warning("blob_classify() : Should probably change from x,y to E,S at some point...");
            % Turn `isroof` and `issurface` into proper masks
            obj.isroof = double(logical(obj.roofs));
            obj.issurface = double(logical(obj.issurface) + logical(obj.isroof));   % Modified roofs are added back in
            obj.isroof(obj.isroof == 0) = nan;
            obj.issurface(obj.issurface == 0) = nan;
            
            % Indicate which roofs are solarPV-viable (in current state)
            obj.isroof_viable = zeros(size(obj.isroof));
            for i = 1:max(obj.roofs,[],'all')
                if obj.info.type(i) ~= "irregular"
                    obj.isroof_viable = obj.isroof_viable + (obj.roofs==i);
                end
            end
            obj.isroof_viable(obj.isroof_viable == 0) = nan;
            
%             % (!!!TEMP!!!) Add a slope to CCIS Slice 13 (main rooftop)
%             obj.Z = obj.Z + (obj.roofs==13).*((0.3/12)*(obj.X-80)+(-5/12)*(obj.Y-80));
        end
        
        function outputArg = plot(obj,arr)
            % plot(mode) creates a plot of the building.
            %   
            %   Parameters:
            %   -----------
            %   arr : double (2D or 3D array)
            %       Will check if it matches `Z`, `map`, or `map_r`. If
            %       not, then:
            %       - 2D : Assumed to be a height map of the building (see 
            %       `heights_gen`) using same pix-to-m scaling as Z
            %       - 3D : Assumed to be a colour image (e.g. 
            %       `imread('building.png')`) using same pix-to-m scaling
            %       as `map`, `map_r`
            disp_marker = false;
            
            if ndims(arr)==2
                % Read output file:
                output_name = "output/"+obj.name+"_heightmap_"+obj.length+"_percent2m.mat";
                load(output_name,'X1','X2','Y1','Y2','xmax','ymax','theta_CCW') % Loads 'X1','X2','Y1','Y2','xmax','ymax'.
                % Rotate coordinates
                map_dummy = ones(ymax,xmax);
                [map_dummy,coord1,coord2] = map_rotate(map_dummy,theta_CCW,[X1,Y1],[X2,Y2]);
                X1 = coord1(1); Y1 = coord1(2);
                X2 = coord2(1); Y2 = coord2(2);
                % Set up axes, in m
                x_plot = 1:size(arr,2);   y_plot = 1:size(arr,1);   % Not using obj.X, etc., in case arr is different
                [X_plot,Y_plot] = meshgrid(x_plot,y_plot);
                Z_plot = arr;        % Assumed to already be in m
                X_plot = (X_plot-1) * obj.pix2m;
                Y_plot = (Y_plot-1) * obj.pix2m;
                % Create plot
                width=1366;height=768;
                set(gcf,'position',[200,200,width,height]);
                surf(X_plot,Y_plot,Z_plot,'edgecolor','none'); set(gca,'ydir','reverse');
                colorbar; 
                shading interp;
                axis image;  view(2);
                xlabel('X (m)'); ylabel('Y (m)'); colormap(jet(256));
                % Z-limits, titles
                if isequaln(arr,obj.dZdX) || ...
                   isequaln(arr,obj.dZdY) || ...
                   isequaln(arr,sqrt(obj.dZdX.^2 + obj.dZdY.^2))
                    % If it's a gradient map
                    zlim([-1 3]);
                    caxis([-1 3]);
                    daspect([1 1 0.3]);
                    if isequaln(arr,obj.dZdX)
                        title(sprintf("Gradient Map of "+obj.name+", \ndZ/dX"));
                    elseif isequaln(arr,obj.dZdY)
                        title(sprintf("Gradient Map of "+obj.name+", \ndZ/dY"));
                    else
                        title(sprintf("Gradient Magnitude Map of "+obj.name));
                    end
                elseif isequaln(arr,obj.Z)
                    % Assumed to be a height map, with same pixel-to-m
                    % conversions as b.Z.
                    title(sprintf("Height Map of "+obj.name+"\n(m)"));
                else
                    % Other?
                    title(sprintf("Something related to the Height Map of "+obj.name+",\nmaybe?"));
                end
                % Add indicators for the 2 selected points
                if disp_marker
                    z_limits = zlim;
                    hold on; plot3([(X1-1)*obj.pix2m, (X1-1)*obj.pix2m], ...
                                   [(Y1-1)*obj.pix2m, (Y1-1)*obj.pix2m], ...
                                   [z_limits(1), z_limits(2)], ...
                                   'LineWidth',2,'Color','red'); hold off;
                    hold on; plot3([(X2-1)*obj.pix2m, (X2-1)*obj.pix2m], ...
                                   [(Y2-1)*obj.pix2m, (Y2-1)*obj.pix2m], ...
                                   [z_limits(1), z_limits(2)], ...
                                   'LineWidth',2,'Color','red'); hold off;
                end
                rotate3d on;
            elseif ndims(arr)==3
                % Assumed to be a colour map, with same pixel-to-m
                % conversions as b.map/b.map_r.
                % Read output file, for percent positions of 2 points:
                xmax_plot = size(arr,2);    ymax_plot = size(arr,1);
                output_name = "output/"+obj.name+"_map_percent2m.mat";
                load(output_name,'map_r','X1_r','X2_r','Y1_r','Y2_r','xmax_r','ymax_r', ...
                                 'map','X1','X2','Y1','Y2','xmax','ymax', ...
                                 'pix2m_x','pix2m_y');
                if isequal(arr,obj.map)
                    X1_p = (X1-1) / (xmax_plot-1);    X2_p = (X2-1) / (xmax_plot-1);
                    Y1_p = (Y1-1) / (ymax_plot-1);    Y2_p = (Y2-1) / (ymax_plot-1);  
                else %if isequal(arr,obj.map_r), or any other case, then it's assumed to be reprojected
                    X1_p = (X1_r-1) / (xmax_plot-1);    X2_p = (X2_r-1) / (xmax_plot-1);
                    Y1_p = (Y1_r-1) / (ymax_plot-1);    Y2_p = (Y2_r-1) / (ymax_plot-1);   
                end
                % Plot
                imshow(uint8(arr), 'XData', (0:(xmax_plot-1))*obj.pix2m_r_x, 'YData', (0:(ymax_plot-1))*obj.pix2m_r_y); 
                axis on;    xlabel('X (m)');	ylabel('Y (m)');    daspect([1 1 1]);
                title("Colour Map of "+obj.name);
                % Add indicators for the 2 selected points
                if disp_marker
                    hold on; scatter(X1_p*(xmax_plot-1)*obj.pix2m_r_x,Y1_p*(ymax_plot-1)*obj.pix2m_r_y,1000,'red','X'); hold off;
                    hold on; scatter(X2_p*(xmax_plot-1)*obj.pix2m_r_x,Y2_p*(ymax_plot-1)*obj.pix2m_r_y,1000,'red','X'); hold off;
                end
            end
        end
    end
end

