function varargout = heights_gen(varargin)
%HEIGHTS_GEN saves & returns a height map for a building `name`, in pixels.
%Will generate a new height map if one has not yet been created.
%
%   Parameters:
%   -----------
%   name : string
%       Name of building. Case-sensitive.
%   length : integer(?)
%       Desired pixel length of longer side of model. Higher numbers mean
%       better resolution, but longer computing times.
%
%   Outputs:
%   --------
%   `Z = heights_gen(name,length)`:
%       `Z` is a 2D matrix representing the height map of the building.
%   `[Z,X,Y] = heights_gen(name,length)`:
%       `X`,`Y` are meshgrids of Z's coords. In metres if measured.
%   `[Z,X,Y,dZdX,dZdY] = heights_gen(name,length)`:
%       `dZdX`,`dZdY` matrices displaying the slopes of the height map.
%   `[Z,X,Y,dZdX,dZdY,p] = heights_gen(name,length)`:
%       `p` is the mesh, i.e. a `struct` containing the faces, vertices,
%       and normals.

    switch nargin
        case 1
            name = varargin{1};
            warning("Building(): No pixel length specified; setting to default of 400 (no rotation).");
            length = 400;
            theta_CCW = 0;
            input_path_meshes = "meshes/";
            output_path_save = "output/";
        case 2
            name = varargin{1};
            length = varargin{2};
            theta_CCW = 0;
            input_path_meshes = "meshes/";
            output_path_save = "output/";
        case 3
            name = varargin{1};
            length = varargin{2};
            theta_CCW = varargin{3};
            input_path_meshes = "meshes/";
            output_path_save = "output/";
        case 5
            name = varargin{1};
            length = varargin{2};
            theta_CCW = varargin{3};
            input_path_meshes = varargin{4};
            output_path_save = varargin{5};
        otherwise
            error('Use `heights_gen(name,length,theta_CCW)`.');
    end
    
%     input_path_meshes = "meshes/";
    fname  = input_path_meshes+name+".stl";
    output_name = output_path_save+name+"_"+length+"pix_"+theta_CCW+"deg.mat";
    
    % Read backup (if any) to check for changes to the STL file
    if isfile(output_name)
        disp('Output found');
        p = nan;
        load(output_name,'p')
        p_backup = p;
    else
        disp('Output not found');
        % No backup found!
        p_backup = nan;
    end
        
    % read STL
    fprintf('Reading STL...\n');
    p = stl_to_pix(fname,length,theta_CCW);
    fprintf('... reading finished!\n');
    faces = p.faces;
    vertices = p.vertices;
    normals = p.normals;
    
    % If "<name>_<length>pix.mat" does not exist, OR "<name>.stl" has been
    % modified since last run, OR "<name>.stl" has not been used before
    if ~isfile(output_name) ...
            | ~isequal(p,p_backup) ...
            | ~isstruct(p_backup) % <-- Redundant?
        % Create backup of new mesh "p", and save it
        % Build triangulation data structure (for Cartesian->Barycentric)
        fprintf('Either the .stl file has been altered or the heightmap has not been generated before. Generating new heightmap...\n');
        TR = triangulation(faces,vertices(:,1:2));     % triangulation(T,P)
                                                       % But ignore z-axis for now.
        % Top-right vertex gives image size
        Z_size = fliplr(ceil( max(vertices(:,1:2)) ));
        rows = Z_size(1);       % max(Y in pixels) = number of rows
        cols = Z_size(2);       % max(X in pixels) = number of cols
        x = 1:cols;
        y = 1:rows;
        % outside of triangulation you have NAN
        Z = nan( [rows,cols] );
        dZdX = nan( [rows,cols] );
        dZdY = nan( [rows,cols] );

        % Begin looping through faces!
%         fprintf('Starting loop:\n');
        for i = 1:size(faces,1)
            face = faces(i,:);
            norm = normals(i,:);

            % get this triangle's bounding box   
            vs = vertices(face,1:2);        % Each row corresponds to one vertex in the face.
            xymin = floor( min(vs(:,1:2) )); % Minimum X,Y values containing triangle. Z-axis (3rd column) ignored.
            xymax = ceil( max(vs(:,1:2)) );  % Maximum X,Y values containing triangle. Z-axis (3rd column) ignored.

            % implicit coordinates for triangle
            [X,Y] = meshgrid( xymin(1):xymax(1), xymin(2):xymax(2) );
            XY = [X(:),Y(:)];            % Each row is a pair of pixel coordinates (in X,Y grid) in bounding box.
                                         % Needs to be in this form to convert each coordinate to barycentric.
            TI = i*ones(size(XY,1),1);   % Current face number (i.e. index of triangle), matching each coord in XY
                                         % to a specific triangle (i) in triangulation TR.

            % get barycentric & check interior
            bary  = TR.cartesianToBarycentric(TI,XY);     % XY, but in barycentric coordinates relative to triangulation TR.
    %         bary2 = cartesianToBarycentric(TR,TI,XY);     % Same as above. This is used in MATLAB's documentation.
            inside = ( all( bary>=0, 2 ) & all( bary<=1, 2 ) ); % 1D boolean array of points in bounding box.
                                                                % Inside triangle (1) if all three vertex coords are >=0 and 
                                                                % <=1. Otherwise 0.
                                                                % ("all()" tests along rows, outputting a column.)
            % set face index inside the image
            X_final = XY(inside,1);
            Y_final = XY(inside,2);
            idxs = sub2ind(size(Z),Y_final,X_final);
            % ^ Syntax: idxs = sub2ind(sz,row,col),
            %   where row, col are (together pair-forming) arrays of subscripts 
            %   ("2D indices") within matrix of size "sz" (in [rows,cols]) that 
            %   you'd like converted to 1D indices.

            % Get final heights!
            p0 = vertices(face(1),1:3);   % x0,y0,z0; one of the vertices. Doesn't matter which.
            Z_newer = p0(3) - 1/norm(3)*(norm(1)*(X_final-p0(1))+norm(2)*(Y_final-p0(2)));  % Heights of this triangle

            % Update heights if "new" heights are higher
            Z_newer_index = Z(idxs)<Z_newer | isnan(Z(idxs));     % =1 where Z<Z_newer, or where Z==nan.
            Z(idxs) = nansum([~Z_newer_index.*Z(idxs), Z_newer_index.*Z_newer],2);
            % Update slopes based on idxs!
            dZdX(idxs) = nansum([~Z_newer_index.*dZdX(idxs), -Z_newer_index.*norm(1)/norm(3)],2);
            dZdY(idxs) = nansum([~Z_newer_index.*dZdY(idxs), -Z_newer_index.*norm(2)/norm(3)],2);

            %Display progress
            if i/5000 == floor(i/5000) | i==size(faces,1)
                disp(""+i/size(faces,1)*100+"% complete");
            end
        end

        % Save output workspace variables to .mat file
        [X,Y] = meshgrid(x,y);
        save(output_name,'Z','dZdX','dZdY','X','Y','p');
    else    % if "<name>_<length>pix.mat" exists, and "<name>.stl" is unchanged
        fprintf('Loading previous results!\n');
        load(output_name);
    end
    
    switch nargout        
        case 3
            varargout{1} = Z;
            varargout{2} = X;
            varargout{3} = Y;
        case 5
            varargout{1} = Z;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = dZdX;
            varargout{5} = dZdY;
        case 6
            varargout{1} = Z;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = dZdX;
            varargout{5} = dZdY;
            varargout{6} = p;
        otherwise
            varargout{1} = Z;
    end
end