function varargout = pix_to_m(varargin)
% function arr = pix_to_m(name,arr)
%PIX_TO_M returns conversion factor between pixels and metres for the input
%   map. Will prompt the user to generate their own if a conversion factor
%   does not already exist.
%
%   Parameters:
%   -----------
%   name : string
%       Name of building. Case-sensitive.
%   arr : double (2D or 3D array)
%       2D : Assumed to be a height map of the building (see `heights_gen`)
%       3D : Assumed to be a colour image (e.g. `imread('building.png')`)
%   theta_CCW (optional) : double
%       Angle that the map has been rotated CCW, in degrees. (Currently
%       used only for height map.)
%   output_path_save (="output/") : string
%       Output "savefile" path.
%   input_path_photos (="images/") : string
%       Input path for photos. Used for getting satellite photo height in
%       metres.
%
%   Returns:
%   --------
%   [pix2m] : double
%       Conversion factor between pixels and metres for x-, y-, and z-axes.
%       Only works when `arr` is a 2D height map (which has
%       equally-proportioned axes).
%   [pix2m_x, pix2m_y] : double
%       Conversion factors between pixels and metres, for x-axis and y-axis
%       respectively.
%       Only works when `arr` is a 3D "colormap" (e.g. satellite photo).

    output_path_save = "output/";
    input_path_photos = "images/";

    switch nargin
        case 2
            name = varargin{1};
            arr = varargin{2};
            theta_CCW = 0;
        case 3
            name = varargin{1};
            arr = varargin{2};
            if isequal(class(varargin{3}),"char")
                theta_CCW = 0;
                output_path_save = varargin{3};
            else
                theta_CCW = varargin{3};
            end
        case 4
            name = varargin{1};
            arr = varargin{2};
            theta_CCW = varargin{3};
            output_path_save = varargin{4};
        case 5
            name = varargin{1};
            arr = varargin{2};
            theta_CCW = varargin{3};
            output_path_save = varargin{4};
            input_path_photos = varargin{5};
        otherwise
            error("Use `pix_to_m(name,arr)` or `pix_to_m(name,arr,theta_CCW)`.");
    end
    switch size(arr,3)
        case 1
            % Assumed to be height map
            xmax = size(arr,2);   ymax = size(arr,1);
            length = max(xmax,ymax);
            output_name_heightmap = output_path_save+name+"_heightmap_"+length+"_percent2m.mat";
            output_name           = output_name_heightmap;
            mode = 'height';
        case 3
            % Assumed to be colour image. Heightmap is used to convert to
            % m, since the colour map isn't perfectly top down.
            files = dir(output_path_save);          % Struct of filenames
            files = join(string( natsortfiles({files.name}) )); % Space-separated filenames in one string, naturally sorted
            expression = name+"_heightmap_\w*_percent2m.mat";
            files_heightmap = regexp(files,expression,"match");
            output_name_heightmap = output_path_save+files_heightmap(end);       % Selects the highest-resolution file.
            output_name           = output_path_save+name+"_map_percent2m.mat";
            mode = 'map';
        otherwise
            output_name_heightmap   = "nope";
            output_name             = "also nope";
            error('Only (N x M x 1) and (N x M x 3) arrays are currently supported.')
    end
    
    if isequal(mode,'height')
        if isfile(output_name_heightmap)
            load(output_name_heightmap,'pix2m','X1','X2','Y1','Y2','xmax','ymax');
            save(output_name_heightmap,'theta_CCW','-append');
            pix2m_x = pix2m;  pix2m_y = pix2m;
        else
            % If output does not exist:
            disp(output_name_heightmap+" not found! Please provide the requested information.");
            Z = arr;
            % Make the plot!
            figure('name',name);
            surf(Z,'edgecolor','none'); set(gca,'ydir','reverse');  colormap(brewermap([],'RdYlBu'));
            colorbar;
            shading interp;
            axis image;
            xlabel('X (pixels)'); ylabel('Y (pixels)'); title(sprintf("Height Map of "+name+", Z\n(pixels)"));
            colormap jet;

            % Set initial window position+size
            width=1366;height=768;
            set(gcf,'position',[200,200,width,height]);

            % Enable input
            view(2); rotate3d on;
            title("Please press ENTER when finished rotating, to begin selecting point A."+newline ...
                +"Press ESC to cancel.");
            waitfor(gcf, 'CurrentCharacter', char(13));

            if ~isempty(findobj('type','figure','name',name))
                title("Select point A.");
                % If figure wasn't closed, input X1,Y1.
                [X1,Y1] = ginput(1);
                set(gcf, 'CurrentCharacter', ' ');
                rotate3d on;
                title("Please press ENTER when finished rotating, to begin selecting point B.");
                waitfor(gcf, 'CurrentCharacter', char(13));
                if ~isempty(findobj('type','figure','name',name))
                    % If figure wasn't closed, input X2,Y2.
                    title("Select point B.");
                    [X2,Y2] = ginput(1);
                    if ~isempty(findobj('type','figure','name',name))
                        % Check if satellite photo height (in m) is
                        % saved in a .mat file. If not, then user will
                        % be prompted to input true distance directly.
                        fname  = input_path_photos+string(name);    % Filename for image AND saved height
                        % Read file
                        if isfile(fname+".mat")
                            % Automatic conversion is possible!
                            warning("pix_to_m() : Automatic conversion to real-world units is possible! Returning bogus `pix2m` for now.");
                            dist_m = nan;       % Bogus distance between points A and B.
                        else
                            % Automatic conversion is not possible, and
                            % needs to be inputted manually.
                            % Next prompt if figure wasn't closed:
                            prompt = "Real horizontal distance between the points, from Google Earth, in metres:";
                            dlgtitle = 'Input';
                            dims = [1 40];
                            answer = inputdlg(prompt,dlgtitle,dims);
                            dist_m       = str2num(answer{1});  % Real distance between points A and B.
                        end
                        dist_pix     = sqrt( (X2-X1)^2 + (Y2-Y1)^2 );   % Distance between selected points on hmap in "pixels"
                        pix2m     = dist_m / dist_pix;
                        pix2m_x   = pix2m;
                        pix2m_y   = pix2m;
                        xmax = size(Z,2);   ymax = size(Z,1);
                        close;
                        % Convert to m; display once more.
                        % (These arrays will start at ZERO, not 1, since they 
                        % don't match indices anymore.)
%                         x = 1:xmax; y = 1:ymax;
%                         [X,Y] = meshgrid(x,y);
%                         Z = (Z-1) * pix2m;
%                         X = (X-1) * pix2m;
%                         Y = (Y-1) * pix2m;
%                         surf(X,Y,Z,'edgecolor','none');
%                         set(gca,'ydir','reverse');  colormap(brewermap([],'RdYlBu'));
%                         colorbar;
%                         shading interp;
%                         axis image;
%                         xlabel('X (m)'); ylabel('Y (m)'); title(sprintf("Height Map of "+name+"\n(m)"));
%                         colormap jet;
%                         rotate3d on;
                    end
                end
            end
            save(output_name,'pix2m','X1','X2','Y1','Y2','xmax','ymax','theta_CCW');
        end
    elseif isequal(mode,'map')
        if ~isfile(output_name_heightmap)
            % If 3D output does not exist:
            error("3D mesh needed to account for camera tilt. However, "+...
                 output_name_heightmap+" not found.");
            heights_gen(name);
            % ^ This didn't come up much in the main algorithm, so it's
            % left as an error.
        end
        % Rename the `<building>_heightmap_percent2m.mat` variables for
        % safekeeping!
        load(output_name_heightmap,'pix2m','X1','Y1','X2','Y2','xmax','ymax');

        if isfile(output_name)
            load(output_name,'map_r','X1_r','X2_r','Y1_r','Y2_r','xmax_r','ymax_r', ...
                             'map','X1','X2','Y1','Y2','xmax','ymax', ...
                             'pix2m_x','pix2m_y','rotated_CCW');
            % ^ (for COLOUR MAP)
        else
            % If output does not exist:
            disp(output_name+" not found! Please provide the requested information.");

            % Make the plot!
            figure('name',name);
            map = arr;
            xmax_c = size(map,2); ymax_c = size(map,1);
            imagesc(map); daspect([1 1 1]); %set(gca, 'ydir', 'normal'); 
            axis image;
            xlabel('X (pixels)'); ylabel('Y (pixels)'); title(sprintf("Colour Map of "+name));

            % Set initial window position+size
            width=1366;height=768;
            set(gcf,'position',[200,200,width,height]);

            % Enable input
            zoom on;
            title("Please press ENTER when finished zooming, to begin selecting point A."+newline ...
                +"Press ESC to cancel.");
            waitfor(gcf, 'CurrentCharacter', char(13));

            if ~isempty(findobj('type','figure','name',name))
                title("Select point A.");
                % If figure wasn't closed, input X1,Y1.
                [X1_c,Y1_c] = ginput(1);
                set(gcf, 'CurrentCharacter', ' ');
                zoom on;
                title("Please press ENTER when finished zooming, to begin selecting point B.");
                waitfor(gcf, 'CurrentCharacter', char(13));
                if ~isempty(findobj('type','figure','name',name))
                    % If figure wasn't closed, input X2,Y2.
                    title("Select point B.");
                    [X2_c,Y2_c] = ginput(1);
                end
                zoom on;
            end
            close;
            
                % BEFORE CONTINUING: Convert height map's pixel units to
                % real-world units!
                
                % Check if satellite photo height (in m) is
                % saved in a .mat file. 
                % If not, then it is assumed that pix2m_h has already been
                % found (through the manual Google Earth method).
                fname  = input_path_photos+string(name);    % Filename for image AND saved height
                % Read satellite photo (which SHOULD BE IDENTICAL to `map`, but
                % may be rotated from it depending on code alterations)
                if isfile(fname+".png")
                    photo = imread(fname+".png");
                elseif isfile(fname+".jpg")
                    photo = imread(fname+".jpg");
                else
                    error("pix2m(cmap) : "+fname+".png and "+fname+".jpg do not exist.");
                    return;
                end
                % Read .mat file
                if isfile(fname+".mat")
                    % Automatic conversion is possible!
                    load(fname+".mat","cmap_height");
                    dist_pix_h  = sqrt( (X2-X1)^2 + (Y2-Y1)^2 );            % Distance between selected points on hmap, in "pixels"
                    dist_pix_c    = sqrt( (X2_c-X1_c)^2 + (Y2_c-Y1_c)^2 );  % Distance between selected points on cmap, in pixels
                    c_per_h = dist_pix_c / dist_pix_h;              % Colourmap pixels (before reprojection) per heightmap pixel
                    pix2m_c = str2num(cmap_height) / size(photo,1); % Metres per colourmap pixel (before reprojection)
                    pix2m = pix2m_c * c_per_h;                      % Heightmap pixels
                    save(output_name_heightmap,'pix2m','X1','Y1','X2','Y2','xmax','ymax');
                else
                    % Automatic conversion is not possible, and
                    % needs to be inputted manually.
                    if isnan(pix2m)
                        error("pix2m(cmap) : '"+name+".mat' does not exist so pix2m (for hmap) should have already been inputted manually... but it wasn't!");
                    end
                    % Do nothing
                end
                
            % Swap variable names around: <X> -> <X_h>,
            % and <X_c> -> <X>.
            X1_h = X1;    Y1_h = Y1;
            X2_h = X2;    Y2_h = Y2;
            xmax_h = xmax;  ymax_h = ymax;
            pix2m_h = pix2m;
            clear('X1','Y1','X2','Y2','xmax','ymax','pix2m');    % Renames all of these variables.

            xmax = xmax_c;  ymax = ymax_c;
            X1 = X1_c;      X2 = X2_c;
            Y1 = Y1_c;      Y2 = Y2_c;
                
            % Create dummy heightmap for reprojecting
            heightmap_dummy = zeros(ymax_h,xmax_h,1);
            % Reproject!
            [map_r, X1_r, X2_r, Y1_r, Y2_r, rotated_CCW] = map_reproject(map,heightmap_dummy,X1,X2,Y1,Y2,X1_h,X2_h,Y1_h,Y2_h);
            % Rotate back to normal!
            [map_r,coord1,coord2] = map_rotate(map_r,-rotated_CCW,[X1_r,Y1_r],[X2_r,Y2_r]);
            X1_r = coord1(1); Y1_r = coord1(2);
            X2_r = coord2(1); Y2_r = coord2(2);
            xmax_r = size(map_r,2);   ymax_r = size(map_r,1);
            % Rotate the dummy heightmap as well!
            [heightmap_dummy] = map_rotate(heightmap_dummy,-rotated_CCW);
            xmax_h = size(heightmap_dummy,2);   ymax_h = size(heightmap_dummy,1);

            % Find pixel-to-metre conversions for x- and y-directions
            dist_x = (xmax_h-1)*pix2m_h;   % Same as map it was projected to
            dist_y = (ymax_h-1)*pix2m_h;
            pix2m_x = dist_x / (xmax_r - 1);
            pix2m_y = dist_y / (ymax_r - 1);

            % Display map once more; make sure it's correct!
%             imshow(map_r, 'XData', (0:(xmax_r-1))*pix2m_x, 'YData', (0:(ymax_r-1))*pix2m_y); axis on;
% %             imshow(map_r); axis on;
%             xlabel('X (m)');    ylabel('Y (m)');
%             daspect([1 1 1]);
%             
%             hold on; scatter(X1_r*pix2m_x,Y1_r*pix2m_y,1000,'red','X'); hold off;
%             hold on; scatter(X2_r*pix2m_x,Y2_r*pix2m_y,1000,'red','X'); hold off;
% %             hold on; scatter(X1_r,Y1_r,1000,'red','X'); hold off;
% %             hold on; scatter(X2_r,Y2_r,1000,'red','X'); hold off;
%             title("Reprojected colour map. If incorrect, delete '"+replace(output_name,'_','\_')+"' and retry!")
            save(output_name,'map_r','X1_r','X2_r','Y1_r','Y2_r','xmax_r','ymax_r', ...
                             'map','X1','X2','Y1','Y2','xmax','ymax', ...
                             'pix2m_x','pix2m_y','rotated_CCW');
        end
        if ~isfile(output_name_heightmap)
            % If 3D output does not exist:
            disp("3D mesh needed to account for camera tilt. However, "+...
                 output_name_heightmap+" not found. Running `pix_to_m` for 3D mesh.");
            heights_gen(name);
        end
    end
    
    % Return results  
    switch nargout        
        case 1
            varargout{1} = pix2m;
        case 2
            varargout{1} = pix2m_x;
            varargout{2} = pix2m_y;
        case 3
            varargout{1} = map_r;
            varargout{2} = pix2m_x;
            varargout{3} = pix2m_y;
        otherwise
            error("Can output `[pix2m]`, `[pix2m_x, pix2m_y]`, or `[map_r,pix2m_x, pix2m_y]`.");
    end
    
end

