classdef Building
    % BUILDING object, where name is specified as input.
    %   Contains info about 3D mesh and 2D topview of building, as well as
    %   conversions from pixel to real-life units.
    %
    % Parameters:
    % -----------
    % name : string
    %   Name of building, matching all filenames.
    %
    % 'pixlength' (default: 1600) : string
    %   Number of pixels comprising length of height map.
    %
    % 'panelSpecs' (default: [65,39]) : double
    %   Panel [length,width], in inches (or specified units using
    %   `symunit`).
    %
    % 'panelTilt' (default: 0) : double
    %   Panel tilt (from horizontal), in degrees, for the purpose of
    %   estimating required rooftop space.
    %
    % 'panelOrientation' (default: "landscape") : string
    %   Panel orientation ("landscape" or "portrait").
    %
    % 'walkway_length' (default: 1.8) : double
    %   Minimum required walkway width between panel and edge, in m (or 
    %   specified units using `symunit`).
    %
    % 'force_run_poly' (default: false) : logical
    %   If enabled, forces polygon approx. to run regardless of whether
    %   saved files exist.
    %
    % 'lat' (default: 53.5461 for Edmonton) : double
    %   Latitude of building's location.
    %
    % 'force_run_shadow' (default: false) : logical
    %   If enabled, forces shadow study to run from scratch regardless of
    %   whether saved file exists.
    %
    % 'shadow_scale' (default: 0.3) : double
    %   Resolution (from 0 to 1) of heightmap used for calculating shadows. 
    %   Computation time increases with shadow_scale^2.
    %
    % 'input_path_meshes' (default: 'meshes/') : string
    %   Path to input building meshes.
    %
    % 'input_path_photos' (default: 'images/') : string
    %   Path to input building satellite photos.
    %
    % 'output_path_save' (default: 'output/') : string
    %   Output path containing "save" folders between runs
    %   (heightmaps, polygons, shadow study, etc).
    %
    % 'output_path' (default: 'output/transfer') : string
    %   Output path containing final tables & rooftop information
    %
    
    properties
        name = 'Null'   % Building name (case-sensitive).
        pix_length;     % Pixel length of 3D map (user-specified).
        mesh;           % 3D mesh info (faces, vertices, normals).
        Z;              % Height map.
        Z_nomask;       % Height map, without mask applied.
        x;              % X-coordinate array, in m.
        y;              % Y-coordinate array, in m.
        X;              % X-coordinate meshgrid.
        Y;              % Y-coordinate meshgrid.
        dZdX;           % Slope in X-direction.
        dZdY;           % Slope in Y-direction.
        m;              % Magnitude of slope.
        pix2m;          % Metres/pixel, for Z-map
        terrain;        % Height map of approximate terrain elevation.
%         issurface;      % Binary map indicating non-ground surfaces.
        isroof;         % Binary map indicating roofs.
        isroof_p;       % Binary map indicating roofs, WITH parapets.
        isroof_viable   % Binary map indicating roofs that are either flat or slanted (NOT curved or irregular).
        roofs;          % Blob map (from `bwlabel()`) of each roof.
        roofs_p;        % Blob map (from `bwlabel()`) of each roof, WITH parapets (that weren't removed as a consequence of edge detection).
        roof_isbad;     % Regions where each roof (value corresponding to `roofs`) has unreliable values.
        info;           % Table containing information of each blob.
        shadows;        % Relative amount of shading on building.
        % Trash pile
        map;            % Colour image of the building.
        map_r;          % Colour image of the building, reprojected to mesh's proportions.
        pix2m_r_x;      % Metres/pixel, for reprojected C-map
        pix2m_r_y;      % Metres/pixel, for reprojected C-map
    end
    
    methods
        function obj = Building(name,varargin)
            % CONSTRUCT: Allow for variable inputs.
            % First, establish defaults.
            u = symunit;                            % Allows for easy unit conversions.
            obj.name = name;
            obj.pix_length = 1600;
            panel_orientation = 'landscape';
            panel_tilt = 0*u.deg;                   % Tilt angle of panel from true horizontal (NOT from rooftop)
            panel_unit = u.in;
            panel_y = 65*panel_unit;                % Min. length of rooftop, including regulatory distance from edges
            panel_x = 39*panel_unit;                % Min. width of rooftop, including regulatory distance from edges
            roof_edge_unit = u.m;
            roof_edgedistance = 1.8*roof_edge_unit; % Regulatory distance from edge for each PV unit
            force_run_poly = true;                  % Run the code and forcibly overwrite previous file?
            shadow_scale = 0.3;                     % 0.3 DEFAULT
            lat = 53.5461;                          % Edmonton's latitude (deg, North)
            force_run_shadow = false;               % Run the code and forcibly overwrite previous file?
            
            % Folder names
            input_path_meshes = "meshes/";
            input_path_photos = "images/";
            output_path_save  = "output/";
            output_path       = "output/transfer/";
            
            % Read inputs
            switch nargin
                case 0
                    warning("Building(): No pixel length specified; setting to default of 1600.");
                case 1
                    obj.pix_length = varargin{1};
                otherwise
                    for i = 1:(nargin-1)
                        if isequal(lower(varargin{i}),'pixlength') || isequal(lower(varargin{i}),'pix_length')
                            obj.pix_length = varargin{i+1};
                        elseif isequal(lower(varargin{i}),'panelspecs')
                            panel_specs = varargin{i+1};
                            [~,sUnits] = separateUnits(panel_specs);
                            if sUnits(1) == 1
                                panel_y = panel_specs(1)*u.in;
                                panel_x = panel_specs(2)*u.in;
                            else
                                panel_y = panel_specs(1);
                                panel_x = panel_specs(2);
                            end
                        elseif isequal(lower(varargin{i}),'paneltilt')
                            panel_tilt = varargin{i+1};
                            [~,sUnits] = separateUnits(panel_tilt);
                            if sUnits(1) ~= 1
                            	panel_tilt = val(rewrite(panel_tilt,u.deg));
                            end
                        elseif isequal(lower(varargin{i}),'panel_orientation')
                            panel_orientation = varargin{i+1};
                        elseif isequal(lower(varargin{i}),'walkway_length')
                            roof_edgedistance = varargin{i+1};
                            [~,sUnits] = separateUnits(roof_edgedistance);
                            if sUnits(1) == 1
                                roof_edgedistance = roof_edgedistance(1)*u.m;
                            end
                        elseif isequal(lower(varargin{i}),'force_run_poly')
                            force_run_poly = varargin{i+1};
                        elseif isequal(lower(varargin{i}),'lat') || isequal(lower(varargin{i}),'latitude')
                            lat = varargin{i+1};
                        elseif isequal(lower(varargin{i}),'force_run_shadow')
                            force_run_shadow = varargin{i+1};
                        elseif isequal(lower(varargin{i}),'shadow_scale') || isequal(lower(varargin{i}),'shadow_resolution') || isequal(lower(varargin{i}),'shadow_res')
                            shadow_scale = varargin{i+1};
                            
                        elseif isequal(lower(varargin{i}),'input_path_meshes')
                            input_path_meshes = convertStringsToChars(varargin{i+1});
                            if input_path_meshes(end) ~= '/'
                                input_path_meshes = [input_path_meshes,'/'];
                            end
                        elseif isequal(lower(varargin{i}),'input_path_photos')
                            input_path_photos = convertStringsToChars(varargin{i+1});
                            if input_path_photos(end) ~= '/'
                                input_path_photos = [input_path_photos,'/'];
                            end
                        elseif isequal(lower(varargin{i}),'output_path_save')
                            output_path_save = convertStringsToChars(varargin{i+1});
                            if output_path_save(end) ~= '/'
                                output_path_save = [output_path_save,'/'];
                            end
                        elseif isequal(lower(varargin{i}),'output_path')
                            output_path = convertStringsToChars(varargin{i+1});
                            if output_path(end) ~= '/'
                                output_path = [output_path,'/'];
                            end
                        end
                    end
            end
            
            % Now, with the building named: Initialize!
            obj.map = map_read(obj.name,input_path_photos,output_path_save);   % Reads satellite photo. This cmap faces North.
            fname_stl_rot = output_path_save+obj.name+"_map_percent2m.mat";
            
            % ~~~~~~~~~~~~~~~ GET HMAP+CMAP ~~~~~~~~~~~~~~~~~~~~
            while true
                if isfile(fname_stl_rot)
                    load(fname_stl_rot,'rotated_CCW');
                    theta_CCW = -rotated_CCW;
                else
                    theta_CCW = 0;
                    rotated_CCW = nan;
                end
                [obj.Z,obj.X,obj.Y,obj.dZdX,obj.dZdY,obj.mesh] = heights_gen(obj.name,obj.pix_length,theta_CCW,input_path_meshes,output_path_save);
                obj.m = sqrt( (obj.dZdX).^2 + (obj.dZdY).^2 );  % Magnitude of the slope
                obj.x = obj.X(1,:); obj.y = obj.Y(:,1);
                % Get rid of "inf" values
                obj.Z = obj.Z .* isfinite(obj.Z);

                % Convert from pixels to metres! (These arrays will start at
                % ZERO, not 1, since they don't match indices anymore.)
                obj.pix2m = pix_to_m(obj.name,obj.Z,theta_CCW,output_path_save); 
                % ^ Returns pix-to-m conversion factor for height map.
                % ^ Also saves inputted `theta_CCW` into '<name>_heightmap_<length>_percent2m.mat'.
                obj.Z = (obj.Z-1)*obj.pix2m;    % Just for consistency, even though heights are all relative to ground anyways.
                obj.X = (obj.X-1)*obj.pix2m;
                obj.Y = (obj.Y-1)*obj.pix2m;
                obj.x = obj.X(1,:); obj.y = obj.Y(:,1);
                [obj.map_r, obj.pix2m_r_x, obj.pix2m_r_y] = pix_to_m(obj.name,obj.map,output_path_save);
                % ^ Reprojects cmap to orientation of heightmap, given
                %   theta_CCW, but re-rotates cmap to face north. The
                %   output file, '<name>_map_percent2m.mat', tracks how the
                %   cmap had to be rotated to match the height map.
                %   (`rotated_CCW`)
                % ^ Returns pix-to-m conversion factors for reprojected cmap.
                
                % Pixel conversion factor
                pix = newUnit('pix',obj.pix2m*u.m);
                if theta_CCW == -rotated_CCW
                    % ^ When height map is rotated to face North, just like the cmap 
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
            
            % Mask out unwanted parts of Z
            [Z_mask,mask_id] = roi_mask_gen(obj.name,obj.pix_length,obj.Z,obj.map_r,output_path_save);
            obj.Z_nomask = obj.Z;
            obj.Z = obj.Z .* Z_mask;
            
            % Solar PV parameters (!!! Maybe set this up in its own
            % function / object? !!!)
            switch panel_orientation
                case 'landscape'
                    panel_length = panel_y;
                    panel_width  = panel_x*cos(panel_tilt);
                case 'portrait'
                    panel_length = panel_y*cos(panel_tilt);
                    panel_width  = panel_x;
                otherwise
                    error('`panel_orientation` must be "portrait" or "landscape".');
            end
            % Specify minimum roof dimensions (panel + distance from edges)
            panel_length = rewrite(panel_length, u.pix);
            panel_width  = rewrite(panel_width, u.pix);
            roof_edgedistance = rewrite(roof_edgedistance, u.pix);
            roof_minlength = 2*roof_edgedistance+panel_length;
            roof_minwidth  = 2*roof_edgedistance+panel_width; 
%             roof_mindiameter = min([val(roof_minlength),val(roof_minwidth)]);
            roof_mindiameter = min([val(panel_length),val(panel_width)]);
            % Total minimum area of solarPV-feasible rooftop?
            panel_size_pix   = ceil(val(panel_length*panel_width));     % Number of pixels in each solar panel.
            roof_minsize_pix = ceil(val(roof_minlength*roof_minwidth)); % Number of pixels in smallest panel-legal rooftop.
            % Note: `val(x)` is just `double(separateUnits(x))`.    
            
            % ~~ Roof Identification ~~
            % Identify any region that isn't the ground, an edge, or tiny
            m_max = 12/12;
%             obj.issurface = find_surfaces(obj.Z,obj.terrain,obj.pix2m, m_max, roof_minsize_pix);
%             % Identify any surface that can be treated as a roof
%             obj.isroof = find_roofs(obj.issurface, roof_mindiameter);  % Surfaces that may fit an SPV (??? is `blob_clean` doing it correctly?)
            [~,obj.isroof] = find_surfaces(obj.Z,obj.terrain,obj.pix2m, m_max, roof_minsize_pix, roof_mindiameter/2);
            obj.roofs = blob_heights_split(obj.isroof,obj.Z);     % Roofs, properly separated from each other
            obj.isroof = logical(obj.roofs);                      % Binary for the above
            obj.roof_isbad = blob_clean_spikes(obj.isroof, obj.Z); % Roofs, with sudden "spikes" in height cleaned

%             % Remove roofs from `surfaces` (re-add after `roof` modifications)
%             obj.issurface = obj.issurface - obj.isroof;
            
            % Classify roofs & split some `irregular` roofs!
            [obj.roofs,obj.roofs_p,obj.roof_isbad,obj.info] = blob_classify(obj.roofs,obj.roof_isbad,obj.Z,obj.dZdX,obj.dZdY,obj.pix2m,roof_mindiameter,1);
%             [obj.roofs,obj.roof_isbad,obj.info] = blob_classify_noslopesplit(obj.roofs,obj.roof_isbad,obj.Z,obj.dZdX,obj.dZdY,obj.pix2m,roof_mindiameter,1);
            warning("blob_classify() : stdev threshold for plane fitting arbitrarily set to 0.09, after working with CAM's wavy rooftop. Need to test some more on sloped surfaces and irregular surfaces!");
            warning("blob_classify() : Should probably change from x,y to E,S at some point...");
            % Turn `isroof`, `isroof_p`, ~and `issurface`~ into proper masks
            obj.isroof = double(logical(obj.roofs));
            obj.isroof(obj.isroof == 0) = nan;
            obj.isroof_p = double(logical(obj.roofs_p));
            obj.isroof_p(obj.isroof_p == 0) = nan;
%             obj.issurface = double(logical(obj.issurface) + logical(obj.isroof));   % Modified roofs are added back in
%             obj.issurface(obj.issurface == 0) = nan;
            
            % Indicate which roofs are solarPV-viable (in current state)
            obj.isroof_viable = zeros(size(obj.isroof));
            for i = 1:max(obj.roofs,[],'all')
                if obj.info.type(i) ~= "irregular"
                    obj.isroof_viable = obj.isroof_viable + (obj.roofs==i);
                end
            end
            obj.isroof_viable(obj.isroof_viable == 0) = nan;
            
            % Add polygon approximation!
            kmeans_mode = "default";    % Generates each polygon many times under default K-means settings, and picks the best one.
%             kmeans_mode = "simple";     % Generates each polygon once, with a higher K-means "Replicates" setting (for consistency).
%                                         %   MUCH faster, but tends to be a bit less accurate.
%             kmeans_mode = "smart";      % Generates each polygon once, but reruns K-means numerous times to enhance accuracy.
            table_polygons = poly_approx(obj.name,obj.info,obj.roofs,obj.roofs_p,mask_id,kmeans_mode,force_run_poly,output_path_save);
            obj.info = [obj.info table_polygons];

            % Add shading analysis!
            obj.shadows = shading_analysis(obj,lat,obj.pix_length,shadow_scale,force_run_shadow,output_path_save);
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
                output_name = "output/"+obj.name+"_heightmap_"+obj.pix_length+"_percent2m.mat";
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
%                 xlabel('X (m)'); ylabel('Y (m)'); colormap(jet(256));
                xlabel('X (m)'); ylabel('Y (m)'); colormap(brewermap([],'RdYlBu'));
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
        
        function outputArg = roofplot(obj,varargin)
            % Plots the inputted roof, along with its final polygon.
            switch (nargin-1)
                case 1
                    roof_id = varargin{1};
                    mustBeNumericOrLogical(roof_id);
                    constrict_option = true;
                case 2
                    roof_id             = varargin{1};
                    constrict_option = varargin{2};
                    mustBeNumericOrLogical(roof_id);
                    mustBeNumericOrLogical(constrict_option);
                otherwise
                    error('Use `roofplot(roof (int), collapse (logical)). The latter is optional.');
            end
            
            % Convert roof to 2D logical, if needed
            if isequal(size(roof_id),[1,1])
                roof   = (obj.roofs==roof_id);
                roof_p = (obj.roofs_p==roof_id);
            else
                error("Building.roofplot() - Input roof must be an integer corresponding to the roof ID in Building.info.");
            end
            % "Constrict" roof to its own boundaries, if needed
            if constrict_option
                [rows, cols] = find(roof_p);
                edge_buffer = 50;
                row1 = min(rows) - edge_buffer; row2 = max(rows) + edge_buffer;
                col1 = min(cols) - edge_buffer; col2 = max(cols) + edge_buffer;
                row1 = max(row1,1);     row2 = min(row2,size(roof,1));
                col1 = max(col1,1);     col2 = min(col2,size(roof,2));
                roof   = roof(row1:row2,col1:col2);      % Cropped array containing the roof
            else
                [rows, cols] = find(roof_p);
                row1 = 1;               row2 = size(roof,1);
                col1 = 1;               col2 = size(roof,2);
                roof   = roof(row1:row2,col1:col2);      % Cropped array containing the roof
            end
            
            
            % Get polygon
            corners = obj.info.polygons{roof_id};
            confidence = obj.info.polygons_confidence(roof_id);
            buffer_x = min(cols) - col1;
            buffer_y = min(rows) - row1;
            corners_x = corners(:,1) + buffer_x;    % Temporary, for drawing polygon with new window size
            corners_y = corners(:,2) + buffer_y;    % Temporary, for drawing polygon with new window size
            
            n_corners = length(corners);
            map_polygon = zeros(size(roof));
            % Check if polygon approximation was successful.
            if sum(isnan(corners))==0
                for j = 1:n_corners
                    jp1 = rem(j+0,n_corners)+1;  % Next corner
                    x_current = corners_x(j);
                    y_current = corners_y(j);
                    % corner_next
                    x_next = corners_x(jp1);
                    y_next = corners_y(jp1);

                    % Draw a line from current corner to next corner!
                    [x_b,y_b] = bresenham(x_current,y_current,x_next,y_next);
                    x_b(1) = [];    y_b(1) = [];    % Ignore first point.
                    line_len = length(x_b);
                    for k = 1:line_len
                        x_line = x_b(k);
                        y_line = y_b(k);
                        map_polygon(y_line,x_line) = map_polygon(y_line,x_line) + 1;
                    end
                end
                map_polygon = logical(map_polygon);
                % Thicken line if plotting window is too "zoomed out"
                if ~constrict_option
                    map_polygon = imdilate(map_polygon,strel('disk',1));
                end

                % Plot!
                imagesc(map_polygon + roof*0.1, 'YData',-buffer_y+1:size(roof,1)-buffer_y+1, 'XData',-buffer_x+1:size(roof,2)-buffer_x+1); axis image;
    %             imagesc(map_polygon + roof*0.1); axis image; colorbar
                width=1366;height=768;
                set(gcf,'position',[200,200,width,height]);
                corners_x = corners(:,1);   % ACTUAL corner locations from table
                corners_y = corners(:,2);   % ACTUAL corner locations from table
                hold on
                for ii = 1:length(corners_x)
                    text(corners_x(ii),corners_y(ii),num2str(ii),'Color','g','fontsize',20)  % labels are for position
                end
                hold off
                title(obj.name + ", roof "+roof_id+" (polygon confidence: "+confidence+")" ...
                      + newline + "(Light blue = rooftop. Yellow = polygon.)")
                xlabel("X (pixels, relative to smallest box containing rooftop)");
                ylabel("Y (pixels, relative to smallest box containing rooftop)");
            else
                % Plot!
                imagesc(roof*0.1, 'YData',-buffer_y:size(roof,1)-buffer_y, 'XData',-buffer_x:size(roof,2)-buffer_x); axis image;
    %             imagesc(map_polygon + roof*0.1); axis image; colorbar
                width=1366;height=768;
                set(gcf,'position',[200,200,width,height]);
                corners_x = corners(:,1);   % ACTUAL corner locations from table
                corners_y = corners(:,2);   % ACTUAL corner locations from table
                text(-buffer_x,-buffer_y+0.03*size(roof,1),"Polygon approximation failed.",'Color','g','fontsize',20)  % labels are for position
                title(obj.name + ", roof "+roof_id+" (polygon confidence: "+confidence+")" ...
                      + newline + "(Yellow = rooftop.)")
                xlabel("X (pixels, relative to smallest box containing rooftop)");
                ylabel("Y (pixels, relative to smallest box containing rooftop)");
            end
        end
    end
end

