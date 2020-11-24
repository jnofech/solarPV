function new_hmap = heights_check(varargin)
%HEIGHTS_CHECK checks whether a new heightmap has been read. If so, then a
%new mask is needed.
%
%   Parameters:
%   -----------
%   name : string
%       Name of building. Case-sensitive.
%   length : integer(?)
%       Desired pixel length of longer side of model. Higher numbers mean
%       better resolution, but longer computing times.

    switch nargin
        case 1
            name = varargin{1};
%             warning("Building(): No pixel length specified; setting to default of 400 (no rotation).");
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
%         disp('Output found');
        p = nan;
        load(output_name,'p')
        p_backup = p;
    else
%         disp('Output not found');
        % No backup found!
        p_backup = nan;
    end
        
    % read STL
%     fprintf('Reading STL...\n');
    p = stl_to_pix(fname,length,theta_CCW);
%     fprintf('... reading finished!\n');
    
    % If "<name>_<length>pix.mat" does not exist, OR "<name>.stl" has been
    % modified since last run, OR "<name>.stl" has not been used before
    if ~isfile(output_name) ...
            || ~isequal(p,p_backup) ...
            || ~isstruct(p_backup) % <-- Redundant?
        % Create backup of new mesh "p", and save it
        % Build triangulation data structure (for Cartesian->Barycentric)
%         fprintf('Either the .stl file has been altered or the heightmap has not been generated before. Generating new heightmap...\n');
        new_hmap = true;
    else    % if "<name>_<length>pix.mat" exists, and "<name>.stl" is unchanged
        new_hmap = false;
    end
end