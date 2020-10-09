function varargout = crop(varargin)
% map_cropped = crop(map_1,map_2,...,map_n,edge_buffer (optional))
% 
% Crops a 2D map (`logical` or `double`) to the smallest region containing 
% foreground, plus an edge buffer around it. If multiple maps are inputted,
% the size of the cropping will be determined by the first map.
%
% Parameters:
% -----------
% map_n : logical or double (2D)
%   Doubles will treat regions that are not "NaN" or zero as foreground.
%   Multiple maps may be inputted at once, provided they are the same size.
%   All maps will be cropped based on the outer edges of the first map.
% edge_buffer (optional) : double
%   Integer number of pixels around the foreground, if possible. Default=50
%
% Returns:
% --------
% map_cropped : logical or double (2D; smaller than `map`)
%   Cropped version of input.

%   % OLD version: only one map at a time        
%     switch (nargin)
%         case 1
%             map = varargin{1};
%             edge_buffer = 50;
%             mustBeNumericOrLogical(map);
%         case 2
%             map = varargin{1};
%             edge_buffer = varargin{2};
%             mustBeNumericOrLogical(map);
%             mustBeNumeric(edge_buffer);
%         otherwise
%             error('');
%     end
%     map_binary = (map~=0) .* isfinite(map);
% 
%     [rows, cols] = find(map_binary);
%     row1 = min(rows) - edge_buffer; row2 = max(rows) + edge_buffer;
%     col1 = min(cols) - edge_buffer; col2 = max(cols) + edge_buffer;
%     row1 = max(row1,1);     row2 = min(row2,size(map,1));
%     col1 = max(col1,1);     col2 = min(col2,size(map,2));
%     map_cropped   = map(row1:row2,col1:col2);      % Cropped array containing the roof

    %   NEW version: As many maps as you please!
    % Check if edge buffer was specified
    if isequal(size(varargin{end}), [1,1])
        n_maps = nargin-1;
        edge_buffer = varargin{nargin};
        mustBeNumeric(edge_buffer);
    else
        n_maps = nargin;
        edge_buffer = 50;
    end
    assert(n_maps>=1);
    maps = varargin(1:n_maps);
    
    % Get size of final map; measure rows and cols to do the cropping at
    map = maps{1};
    mapsize = size(map);
    map_binary = (map~=0) .* isfinite(map);
    [rows, cols] = find(map_binary);
    row1 = min(rows) - edge_buffer; row2 = max(rows) + edge_buffer;
    col1 = min(cols) - edge_buffer; col2 = max(cols) + edge_buffer;
    row1 = max(row1,1);     row2 = min(row2,size(map,1));
    col1 = max(col1,1);     col2 = min(col2,size(map,2));
    
    % Loop through maps and crop 'em!
    varargout = cell(1,n_maps);     % Preallocate output cell array
    for i = 1:n_maps
        map = maps{i};
        if ~isequal(size(map),mapsize)
            error("All inputted maps must be of the same size! Map #"+i+" is wrong size.");
        end

        map_cropped   = map(row1:row2,col1:col2);      % Cropped array containing the roof
        varargout{i} = map_cropped;
    end
end