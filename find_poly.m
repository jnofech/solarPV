function [corners,score] = find_poly(varargin)
%FIND_POLY() approximates a polygon shape for a rooftop
%   using the Hough transform to detect straight lines,
%   and K-means clustering to simplify any lines that
%   are overlapping.
% 
%   Parameters:
%   -----------
%   roof : double (2D array)
%       Binary image indicating roof.
%   roof_p (OPTIONAL) : double (2D array)
%       Binary image indicating roof, with parapet.
%       Determines outer bounds of the image.
%
%   Returns:
%   --------
%   corners : double (2D array)
%       An n-by-2 array (where n is the number of
%       "vertices" on the polygon), in order of
%       connectivity.
%   score : double
%       Relative confidence in corners, on scale from 0 to 1.
    switch nargin
        case 1
            roof = varargin{1};
            roof_p = roof;
            kmeans_mode = "default";
        case 2
            roof = varargin{1};
            if ismember(class(varargin{2}),['string','char'])
                roof_p = roof;
                kmeans_mode = varargin{2};
            else
                roof_p = varargin{2};
                kmeans_mode = "default";
            end
        case 3
            roof = varargin{1};
            roof_p = varargin{2};
            if ismember(class(varargin{3}),['string','char'])
                kmeans_mode = varargin{3};
            else
                error("`kmeans_mode` (varargin{3}) inputted incorrectly");
            end            
        otherwise
            error("find_poly() : Can only input `roof`, `roof_parapet`, and `kmeans_mode`.");
    end
    % Cut out useless parts of image
    slice_full = roof_p;
    [rows, cols] = find(slice_full);
    row1 = min(rows) - 50; row2 = max(rows) + 50;
    col1 = min(cols) - 50; col2 = max(cols) + 50;
    row1 = max(row1,1);    row2 = min(row2,size(roof,1));
    col1 = max(col1,1);    col2 = min(col2,size(roof,2));
    roof_small   = roof(row1:row2,col1:col2);      % Cropped array containing each blob, EXCLUDING "bad" regions

    % Get outer edges
    roof_f = imfill(roof_small,'holes');    % Remove obstacles within roof.
    edges = bwperim(roof_f);
    
    % Establish bogus "placeholder" that will indicate a failed run
    x_placeholder = [1, 1, size(edges,2), size(edges,2)];
    y_placeholder = [1, size(edges,1), size(edges,1), 1];
    corners_placeholder = [x_placeholder',y_placeholder'];
    
    % Find accuracy of polygon!
    % Polygon will be multiplied by a thickened map of "edges", with each
    %   polygon pixel scoring up to 1 point for max accuracy.
    edges_minscore = 0.5;    % Minimum score, at max radius around edges.
    edges_maxscore = 1;      % Maximum score, within min radius around edges.
    edges_voidscore = -1;    % Negative score, in the "void" wherever edges certainly aren't.
    edges_minradius  = min(round(length(edges)/190),3);
    edges_maxradius  = max(edges_minradius+1,round(length(edges)/55));
    edges_voidradius = round(edges_maxradius/2); % Distance between "max radius (minscore)" and "void (negative score)".
    edges_glow = zeros(size(edges));
    for radius=edges_minradius:edges_maxradius
        edges_glow = edges_glow + imdilate(edges,strel('disk',radius));
    end
    edges_glow(edges_glow==0) = NaN;
    edges_min = nanmin(edges_glow,[],'all');
    edges_max = nanmax(edges_glow,[],'all');
    edges_glow = (edges_glow - edges_min)/(edges_max - edges_min) * (edges_maxscore - edges_minscore) + edges_minscore;
    edges_glow(isnan(edges_glow)) = 0;
    edges_void = ~imdilate(logical(edges_glow),strel('disk',edges_voidradius));
    edges_glow(edges_void) = edges_voidscore;
    
    % At this point, `edges_glow` is a 2D map ranging from edge_minscore to
    % edge_maxscore based on proximity to `edges`, and zero elsewhere.
    
    % The code will attempt a polygon approximation `n_attempts` times.
    % If a final "score" of over 95% is reached, the code will
    % automatically use that as the final result!
    % Otherwise, it will simply select the best result out of those
    % `n_attempts` runs.
    
    if kmeans_mode=="default"
        n_attempts      = 15;
    elseif kmeans_mode=="simple"
        n_attempts      = 1;
    elseif kmeans_mode=="smart"
        n_attempts      = 1;
    else
        warning("find_poly.m : `kmeans_mode` not set to 'default' or 'simple'! Set to 'default' by default.");
        kmeans_mode = "default";
        n_attempts = 15;
    end
    corner_outputs  = cell(n_attempts,2);   % Will hold the 2D output `corners` arrays.
    best_attempt = NaN;     % Will indicate the highest score.
    failed_counter = 0;     % Will indicate the number of failed attempts (i.e. loop counters maxed out)
    for attempt = 1:n_attempts
        corners = corners_placeholder;
        loop_counter = 0;
        if attempt==1
            % First, try simple convex corner detection+sorting
%             disp("Attempt convex corner detection first!");
            corners_switched = pgonCorners(roof_small,10);
            corners = corners_switched;
            corners(:,1) = corners_switched(:,2);
            corners(:,2) = corners_switched(:,1);
        else
            % Hough Line -> KMeans polygon approximation
            while (isequal(corners,corners_placeholder)|| length(corners)==2) && (loop_counter<=30) 
                % The Hough algorithm will crash if the roof is too tiny.
                % This will result in `corners = [NaN,NaN]`.
                corners = find_poly_attempt(roof,roof_p,kmeans_mode);
                loop_counter = loop_counter + 1;    % Loops will increase if runs fail partway through (i.e. `corners` is set to bogus "placeholder" value during attempt).
            end
            if loop_counter >= 30
                failed_counter = failed_counter + 1;
            end
        end

        if isequaln(corners,[NaN,NaN])
            % If `corners` is set to [NaN,NaN] during
            % `find_poly_attempt()`, it's because the roof is too small for
            % Hough line detection.
            corners = corners_placeholder;
        end
%         disp("loop_counter = "+loop_counter);

        % Draw polygon traced by vertices!
        map_polygon = zeros(size(edges));
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

        % Calculate "score" of similarity to actual edge map!
        score = sum(map_polygon .* edges_glow,'all')/sum(edges,'all');
%             imagesc(edges_glow); axis image; colorbar
%             hold on
%             plot( corners(:,2),corners(:,1),'yo','MarkerFaceColor','r',...
%                                             'MarkerSize',11,'LineWidth',2);
%             hold off
        if length(corners)==2
            % "Straight line" polygons receive a score of zero.
            score = 0;
        end
        
        % Save the output and score
        corner_outputs{attempt,1} = corners;
        corner_outputs{attempt,2} = score;

%                 figure; imagesc(map_polygon + 0.2*edges); axis image; colorbar; title("Score: "+score);
% %                 figure; imagesc(2*map_polygon + edges_glow); axis image; colorbar; title("Score: "+score);
%                 width=1366;height=768;
%                 set(gcf,'position',[200,200,width,height]);
%                 x = corners(:,1);
%                 y = corners(:,2);
%                 hold on
%                 for ii = 1:length(x)
%                     text(x(ii),y(ii),num2str(ii),'Color','g','FontSize',18)
%                 end
%                 hold off
        
        % Any score over 0.95 is good enough!
        if score >= 0.95
            best_attempt = attempt;
            break
        end
        % If simple convex corner detection was needed: Cancel early.
        if failed_counter >= 2
            disp("find_poly() : Too many failed attempts. Ending.");
            break
        end
    end
    if isnan(best_attempt)
        % If a score of >0.95 has not been found yet
        scores = [corner_outputs{:,2}];
        best_score = max(scores);
        best_attempt = find(scores==best_score);
    end
    corners = corner_outputs{best_attempt,1};
    score   = corner_outputs{best_attempt,2};
    
    % Draw best polygon!
    map_polygon = zeros(size(edges));
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
    

%             % Smooth out edges
%             edges_s = imclose(edges,strel('disk',10));
%             edges_s = imfill(edges_s,'holes');    % Remove obstacles within roof.
%             edges_s = bwperim(edges_s);
%             figure; imagesc(map_polygon + 0.2*edges + 0.1*edges_s); axis image; colorbar; title("Score: "+score);
%             width=1366;height=768;
%             set(gcf,'position',[200,200,width,height]);
%             x = corners(:,1);
%             y = corners(:,2);
%             hold on
%             for ii = 1:length(x)
%                 text(x(ii),y(ii),num2str(ii),'Color','g','FontSize',18)
%             end
%             hold off

    % If polygon approximation completely failed
    if isequal(corners,corners_placeholder) || length(corners)==2
        disp("find_poly() : Polygon approximation failed.")
        corners = [NaN, NaN];
        score = 0;
    end
            
    % Return to normal "building" coordinates
    corners(:,1) = corners(:,1) + col1;
    corners(:,2) = corners(:,2) + row1;
end

function corners = find_poly_attempt(varargin)
%FIND_POLY_ATTEMPT() approximates a rooftop polygon
%   using the Hough transform to detect straight lines.
%   Multiple overlapping straight lines are 
%   Since K-means has somewhat randomized behaviour,
%   this function may not have satisfactory results.
% 
%   Parameters:
%   -----------
%   roof : double (2D array)
%       Binary image indicating roof (1).
%
%   Returns:
%   --------
%   x_output : double (n-by-1)
%       x-coordinates for the n vertices.
%   y-output : double (n-by-1)
%       y-coordinates for the n vertices.
%   corners : double (2D array)
%       An n-by-2 array (where n is the number of
%       "vertices" on the polygon), in order of
%       connectivity.
    switch nargin
        case 1
            roof = varargin{1};
            roof_p = roof;
            kmeans_mode = "default";
        case 2
            roof = varargin{1};
            if ismember(class(varargin{2}),['string','char'])
                roof_p = roof;
                kmeans_mode = varargin{2};
            else
                roof_p = varargin{2};
                kmeans_mode = "default";
            end
        case 3
            roof = varargin{1};
            roof_p = varargin{2};
            if ismember(class(varargin{3}),['string','char'])
                kmeans_mode = varargin{3};
            else
                error("`kmeans_mode` (varargin{3}) inputted incorrectly");
            end            
        otherwise
            error("find_poly_attempt() : Can only input `roof`, `roof_parapet`, and `kmeans_mode`.");
    end
    if kmeans_mode=="default"
        kmeans_replicates = 1;
%         disp("DEFAULT");
    elseif kmeans_mode=="simple"
        kmeans_replicates = 500;
        disp("SIMPLE");
    elseif kmeans_mode=="smart"
        kmeans_replicates = 1;
        disp("SMART");
    else
        warning("find_poly_attempt() : kmeans_mode not set to 'default' or 'simple'. Setting to 'default' by default");
        kmeans_mode = "default";
        kmeans_replicates = 1;
    end
    
    % Cut out useless parts of image (Is this necessary?)
    slice_full = roof;
    pix_length = length(roof);
    [rows, cols] = find(roof_p);
    row1 = min(rows) - 50; row2 = max(rows) + 50;
    col1 = min(cols) - 50; col2 = max(cols) + 50;
    row1 = max(row1,1);    row2 = min(row2,size(roof,1));
    col1 = max(col1,1);    col2 = min(col2,size(roof,2));
    roof   = slice_full(row1:row2,col1:col2);      % Cropped array containing each blob, EXCLUDING "bad" regions
    edges_length = length(roof);

    % Get outer edges
    roof_f = imfill(roof,'holes');    % Remove obstacles within roof.
    edges = bwperim(roof_f);

    % Smooth out edges
    edge_smooth_radius = round(pix_length / 150);  % Relatively constant scale, regardless of resolution
%     edges_fat = imdilate(edges,strel('disk',edge_smooth_radius));
%     edges_s = imerode(edges_fat,strel('disk',edge_smooth_radius));
    
    edges_s = imclose(edges,strel('disk',edge_smooth_radius));
    
    edges_s = imfill(edges_s,'holes');    % Remove obstacles within roof.
    edges_s = bwperim(edges_s);
    edges = edges_s;
    
        % Find accuracy of a polygon!
        % Polygon will be multiplied by a thickened map of "edges", with each
        %   polygon pixel scoring up to 1 point for max accuracy.
        edges_minscore = 0.5;    % Minimum score, at max radius around edges.
        edges_maxscore = 1;      % Maximum score, within min radius around edges.
        edges_voidscore_outer = 0;     % Negative score, in the "void" wherever edges certainly aren't.
        edges_voidscore_inner = -1;    % Negative score, in the "void" wherever edges certainly aren't.
        edges_minradius  = min(round(length(edges)/190),3);
        edges_maxradius  = max(edges_minradius+1,round(length(edges)/55));
        edges_voidradius = round(edges_maxradius/2); % Distance between "max radius (minscore)" and "void (negative score)".
        edges_glow = zeros(size(edges));
        for radius=edges_minradius:edges_maxradius
            edges_glow = edges_glow + imdilate(edges,strel('disk',radius));
        end
        edges_glow(edges_glow==0) = NaN;
        edges_min = nanmin(edges_glow,[],'all');
        edges_max = nanmax(edges_glow,[],'all');
        edges_glow = (edges_glow - edges_min)/(edges_max - edges_min) * (edges_maxscore - edges_minscore) + edges_minscore;
        edges_glow(isnan(edges_glow)) = 0;
        edges_void = ~imdilate(logical(edges_glow),strel('disk',edges_voidradius));
        edges_void_inner = imclearborder(edges_void);
        edges_void_outer = ~imfill(~edges_void,'holes');
        edges_glow(edges_void_outer) = edges_voidscore_outer;
        edges_glow(edges_void_inner) = edges_voidscore_inner;
        % ^ At this point, `edges_glow` is a 2D map ranging from edge_minscore to
        % edge_maxscore based on proximity to `edges`, and zero elsewhere.

    % Approximate how "long" (as opposed to round/blocky) a rooftop is.
    % i.e. Find the "center of mass" of rooftop minus obstacles, and then
    % find distance between CoM and furthest point on roof edge. Make circle of
    % that radius, and compare that circle's area with the area of the rooftop.
    stats = regionprops(roof);
    centroid = stats.Centroid;
    x_cen = centroid(1);
    y_cen = centroid(2);
    idx = find(edges_s);
    [y_edge,x_edge] = ind2sub(size(edges),idx);      % x- and y-coordinates of edge.

    dist2 = max((y_edge-y_cen).^2 + (x_edge-x_cen).^2);
    idx = find((y_edge-y_cen).^2 + (x_edge-x_cen).^2 == dist2);
    idx = idx(1);           % In case of perfectly equidistant points from centroid
    x_max = x_edge(idx);
    y_max = y_edge(idx);

    r_circ = sqrt((x_max - x_cen)^2 + (y_max - y_cen)^2);
    A_circ = pi*r_circ^2;
    A_roof = sum(roof,'all');           % Roofs with less area compared to A_circ are assumed to be "narrower".


    % Get Hough lines!
    n_lines = 50;               % MAXIMUM number of lines detected.
    linethreshold_pix_correction = (edges_length / 50)*(1 - A_roof / A_circ); % Corrects length threshold for "narrower" structures.
    linethreshold_pix = (edges_length / 13) - linethreshold_pix_correction;   % Number of pixels along the line in order to register as a line.
    
    % The Hough algorithm will crash if the roof is too tiny.
    try
        [H,T,R] = hough(edges_s,'RhoResolution',0.5,'Theta',-90:0.25:89.5);
        %     imagesc(H); axis image; colorbar;
        P  = houghpeaks(H,n_lines,'threshold',linethreshold_pix);           % P (`n_lines`-by-2) contains coords of peaks in Hough (rho,theta) coords.
        h_lines = houghlines(edges_s,T,R,P,'FillGap',5,'MinLength',7);      % Extracts actual lines from P. Each (rho,theta) pair may compose several line segments.
    catch
        warning("Hough line detection crashed due to roof being too small. Running simple convex corner detection.");
        corners = pgonCorners(imfill(edges,'holes'),4);
%         corners = [NaN,NaN];
        return;
    end
    theta = [h_lines.theta]';
    rho = [h_lines.rho]';
    
    % Add perfectly vertical/horizontal lines that may not have
    % been detected.
    %   (Loop through each row (horizontal) and then column (vertical) of `edges`.
    %   If the sum of all elements in that row/col exceeds linethreshold_pix,
    %   then [rho,theta] = [x,0] (vertical) or [-y,-90] (horizontal)!)
    % Check for horizontal lines
    for y = 1:(size(edges_s,1)-1)
        h_sum = sum(edges_s([y,y+1],:),'all');
        if h_sum >= linethreshold_pix
            theta = [theta; -90];
            rho = [rho; -y];
        end
    end
    % Check for vertical lines
    for x = 1:(size(edges_s,2)-1)
        v_sum = sum(edges_s(:,[x,x+1]),'all');
        if v_sum >= linethreshold_pix
            theta = [theta; 0];
            rho = [rho; x];
        end
    end
    lines = [rho,theta];
    
    % K-means: Remove identical lines
    k_clusters = size(unique(lines,'rows'),1);         % Aka "k".
    try
        [~,kmeans_C] = kmeans(lines,k_clusters,'Replicates',kmeans_replicates);% Clusters line segments based on (rho,theta) alone, with random starting points.
    catch
        warning("K-means clustering crashed.");
        corners = [NaN,NaN];
        return;
    end
    % ^ `lines` is an n-by-2 array of [rho,theta] pairs.
    % ^ kmeans_idx is an n-by-1 vector of cluster indices of each observation(???).
    % ^ kmeans_C is a k-by-2 array, indicating the `k` centroid locations in rho-theta space
    lines = kmeans_C;

    % Draw lines
    overlap_count = 0;
    overlap_threshold = 30;     % Number of pixels that a new line can overlap at once before it's considered a duplicate.
    line_maps       = cell(length(lines),1);   % Will display each line on a 2D map.
    line_maps_wide  = cell(length(lines),1);   % Will display thickened versions of each line on a 2D map.
    for i = 1:size(lines,1)    
        [rho_loop]   = lines(i,1);
        [theta_loop] = lines(i,2);
        A = cosd(theta_loop);
        B = sind(theta_loop);
        x0 = A * rho_loop;
        y0 = B * rho_loop;
        x1 = int16(x0 + 10000 * (-B));
        y1 = int16(y0 + 10000 * A);
        x2 = int16(x0 - 10000 * (-B));
        y2 = int16(y0 - 10000 * A);
        % Make sure coords are within image
        pixels = createLineIterator([x1, y1], [x2, y2], edges_s);
        pixels = int16(pixels);
        pixel_len = length(pixels);
        % Save Bresenham line
        line_maps{i}   = zeros(size(edges));
        line_maps_wide{i}   = zeros(size(edges));
        for j = 1:pixel_len
            x = pixels(j,1);
            y = pixels(j,2);
            line_maps{i}(y,x) = line_maps{i}(y,x) + 1;
        end
    end
    % Loop through lines; detect "overlap"
    line_widen_radius = 1;     % Radius of structuring element that widens each line.
    strel_widen = strel('disk',line_widen_radius);
    map_lines = zeros(size(edges));
    for i = 1:size(lines,1)
        line_map = line_maps{i};
        line_maps_wide{i} = imdilate(line_map,strel_widen);
        line_map_wide = line_maps_wide{i};

        map_lines_backup = map_lines;   % Before the new line is added.
        map_lines = map_lines + line_map_wide;

        n_intersections_score      = sum(map_lines(map_lines > 1)-1) / (2*line_widen_radius + 1)^2; % 1 point for every intersection; multiple overlaps count for more points.
        n_intersections_score_prev = sum(map_lines_backup(map_lines_backup > 1)-1) / (2*line_widen_radius + 1)^2;
        if abs(n_intersections_score - n_intersections_score_prev) > overlap_threshold
            overlap_count = overlap_count + 1;
        end
    end

    % Run K-means again, to get final lines!
    if kmeans_mode == "smart"
        % Number of lines VARIES. For each number of lines, run K-means,
        % draw lines, and calculate a final "score". Highest score gets
        % used as actual polygon!
        k_clusters_base = length(unique(lines,'rows')) - overlap_count;
        k_range_width = 3;
        k_clusters_min = max( [3,k_clusters_base-k_range_width,ceil(k_clusters_base*2/3)] );    % Typically k_clusters_base-k_range_width
        k_clusters_max = min( [k_clusters_base+k_range_width,ceil(k_clusters_base*4/3)] );      % Typically k_clusters_base+k_range_width
        k_values = k_clusters_min:k_clusters_max;   % The k-values to be looped through
        k_scores = zeros(size(k_values));   % Maximum score of each `k` over `n_attempts` attempts
        k_lines = cell(size(k_values));     % (will be) Highest scoring set of lines for each `k`
        for k = 1:length(k_values)
            % Number of attempts
            n_attempts = 15;
            % Run k-means `n_attempts` times. Pick best score!
            for attempt = 1:n_attempts
                [~,kmeans_C] = kmeans(lines,k_values(k));           % Clusters line segments based on (rho,theta) alone, with random starting points.
                % Find score
                line_maps_temp = cell(length(kmeans_C),1);   % Will display each line on a 2D map.
                for i = 1:size(kmeans_C,1)    
                    [rho_loop]   = kmeans_C(i,1);
                    [theta_loop] = kmeans_C(i,2);
                    A = cosd(theta_loop);
                    B = sind(theta_loop);
                    x0 = A * rho_loop;
                    y0 = B * rho_loop;
                    x1 = int16(x0 + 10000 * (-B));
                    y1 = int16(y0 + 10000 * A);
                    x2 = int16(x0 - 10000 * (-B));
                    y2 = int16(y0 - 10000 * A);
                    % Make sure coords are within image
                    pixels = createLineIterator([x1, y1], [x2, y2], edges_s);
                    pixels = int16(pixels);
                    pixel_len = length(pixels);
                    % Save Bresenham line
                    line_maps_temp{i}   = zeros(size(edges));
                    for j = 1:pixel_len
                        x = pixels(j,1);
                        y = pixels(j,2);
                        line_maps_temp{i}(y,x) = line_maps_temp{i}(y,x) + 1;
                    end
                    % Calculate "score" of similarity to actual edge map!
                    score = sum(line_maps_temp{i} .* edges_glow,'all')/sum(edges,'all');
                    
                    % If it's the highest score for this `k`: Save it!
                    if score > k_scores(k)
                        k_lines{k} = kmeans_C;
                        k_scores(k) = score;
                    end
                end
            end
        end
        % Find the highest-scoring K-value with the lowest number of lines.
        % And use its lines!
        best_idx = find(k_scores==max(k_scores));
        best_idx = best_idx(1); % Just the first.
        lines = k_lines{best_idx};
    else
        % Number of lines = k_previous - overlap_count.
        k_clusters = length(unique(lines,'rows')) - overlap_count;
    %     [~,kmeans_C] = kmeans(lines,k_clusters);           % Clusters line segments based on (rho,theta) alone.
        [~,kmeans_C] = kmeans(lines,k_clusters,'Replicates',kmeans_replicates);% Clusters line segments based on (rho,theta) alone.
        lines = kmeans_C;
    end

    % Find corner locations
    corners_map         = poly_findcorners(roof,lines,edges,edges_s);
    % Sort the corners!
%     [x_output,y_output] = poly_sortcorners(corners_map,edges);
    [x_output,y_output] = poly_sortcorners(corners_map,edges_s);
    corners = [x_output',y_output'];
    
    % Remove any unneeded corners (i.e. the ones lying on a straight line)!
%     corners = poly_cleancorners(corners,edges);
    corners = poly_cleancorners(corners,edges_s);
    
    % Untangle any tangled lines!
%     corners = poly_untangle(corners,edges,1);
    corners = poly_untangle(corners,edges_s,1);
    
    % Remove any unneeded corners (again), in case untangling reintroduced some!
%     corners = poly_cleancorners(corners,edges);
    corners = poly_cleancorners(corners,edges_s);
end
function corners_list = corners_listbetween(p1,p2,n_corners)
% Lists the corners between p1 (earlier) and p2 (later)
% in a list of `n_corners` corners.
% e.g. if p1 = 9 and p2 = 2 in a 10-corner list, then
%      the result will be [9,10,1,2].
    if p1 > p2
        p2 = p2 + n_corners;
    end
    corners_list = p1:1:p2;
    corners_list(corners_list>n_corners) = corners_list(corners_list>n_corners)-n_corners;
end
function d = corner_distance(p1,p2,n_corners)
% Finds the difference (in counting) between two corners
% in a list of `n_corners` corners.
    d = min(abs(p1-p2), abs(min(p1,p2)+n_corners - max(p2)));
end

function corners_map = poly_findcorners(roof,lines,edges,edges_s)
%POLY_FINDCORNERS() finds the locations of vertices that act as "corners"
% for the n-sided roof polygon.
%   Output:
%   -------
%   corners_map : logical (2D)
%       2D map indicating n-gon corners.
    % Initialize arrays
    map_lines    = zeros(size(roof));       % Will show all lines drawn at once.
    line_pixels  = cell(length(lines),1);   % Will hold the n x 3 "pixels" array.
    line_maps    = cell(length(lines),1);   % Will display each line on a 2D map.
    line_corners = cell(length(lines),1);   % Will hold the 2D mask indicating final corners.

    % Find pixels for each line
    for i = 1:size(lines,1)
        [rho_loop]   = lines(i,1);
        [theta_loop] = lines(i,2);
        A = cosd(theta_loop);
        B = sind(theta_loop);
        x0 = A * rho_loop;
        y0 = B * rho_loop;
        x1 = int16(x0 + 10000 * (-B));
        y1 = int16(y0 + 10000 * A);
        x2 = int16(x0 - 10000 * (-B));
        y2 = int16(y0 - 10000 * A);
        % Create Hough line
        pixels = createLineIterator([x1, y1], [x2, y2], edges_s);
        pixel_len = length(pixels);
        % Save Hough line
        pixels = int16(pixels);
        line_pixels{i} = pixels;
        line_maps{i}   = zeros(size(edges));

        for j = 1:pixel_len
            x = pixels(j,1);
            y = pixels(j,2);
            % Displaying where lines are
            map_lines(y,x)    = map_lines(y,x) + 1;
            line_maps{i}(y,x) = line_maps{i}(y,x) + 1;
        end
    end

    % Create "mask" of intersections
    mask_intersections = map_lines >=2;
    % Create final "mask" of corners (most of which will be intersections).
    corners_map = zeros(size(edges));

    % Set up frequently used structuring elements
    strel_overlap = strel('disk',3);    % Used in detecting `edges` near line segments.
    strel_1pix    = strel('square',3);  % Used in reaching beyond line segment by 1 pixel.
    corner_tolerance = 3;
    strel_corner  = strel('disk',corner_tolerance); % Used in detecting corners that intersections failed to outline.

    % Loop through lines; find "corners"
    for i = 1:size(lines,1)
        line_map = line_maps{i};
        line_corners{i} = zeros(size(edges));
        
    %                 % Uncomment and tweak for debugging plots.
    %                 figure; imagesc(edges + edges_s*0.3 + 3*line_map); axis image; colorbar; title("line "+i); 
    %                 hold on; plot(x_finalcorners,y_finalcorners,'o','color','white'); hold off;
    %                 width=1366;height=768;
    %                 set(gcf,'position',[200,200,width,height]);

        % Separate "segments" with intersection mask. Loop through them.
        segments = bwlabel(line_map .* (~mask_intersections),8);
        n_segments = max(segments,[],'all');
        good_segments_idx = NaN(n_segments,1);
        for j = 1:n_segments
            segment = (segments==j);
%             segment_fat = imdilate(segment,strel_overlap);    % Segment and nearby pixels
%             overlap_fraction = abs(sum(segment_fat .* edges,'all') / sum(segment,'all'));
            edges_checkoverlap = imdilate(edges,strel_overlap);
            overlap_fraction = sum(segment .* edges_checkoverlap,'all') / sum(segment,'all');
            
            % Only consider segments that arise from intersections.
            idx1 = find((imdilate(segment,strel_1pix) - segment + line_map)==2);
            if ~isempty(idx1) 
                % If overlap between "fat segment" and "smoothed edge" is within
                % ~15% of segment's length:
                if overlap_fraction > 0.85
                    % Segment is "good". Mark both endpoints as "corners".
                    good_segments_idx(j) = j;
                    line_corners{i} = line_corners{i} + ((imdilate(segment,strel_1pix) - segment + line_map)==2);
                elseif (overlap_fraction > 0.4)
                    % Segment is assumed to be "good" at one end, but not at the
                    % other.
                    good_segments_idx(j) = j;

                    % Find corner closest to roof edge.
                        idx = find((imdilate(segment,strel_1pix) - segment + line_map)==2);
                        [y_corners,x_corners] = ind2sub(size(edges),idx);
                        idx = find(edges);
                        [y_edge,x_edge] = ind2sub(size(edges),idx);

                        % Loop through corners; find closest distance to any point on edge.
                        corner_edgedistances = zeros(length(y_corners),1);
                        for k = 1:length(y_corners)
                            y_corner = y_corners(k);
                            x_corner = x_corners(k);
                            corner_edgedistances(k) = min(sqrt((y_edge - y_corner).^2 + (x_edge - x_corner).^2));
                        end
                    [~,i_closer] = min(corner_edgedistances);
                    x_corner1 = x_corners(i_closer);
                    y_corner1 = y_corners(i_closer);

                    % Find furthest point on line from edge!
                    % Dilate "edges". The 2nd corner will be `corner_tolerance`
                    % pixels closer than the furthest point from corner1 of
                    % (edges_fat .* segment).
                    idx = find(imdilate(edges,strel_corner) .* segment);
                    [y_line,x_line] = ind2sub(size(edges),idx);
%                         disp('AAAA');
%                         y_line
%                         y_corner1
%                         x_line
%                         x_corner1
                    corner1_distances = sqrt((y_line - y_corner1).^2 + (x_line - x_corner1).^2);
                    [~,idx_distances] = sort(corner1_distances);
                    x_corner2 = x_line(idx_distances(max(1,end-corner_tolerance)));
                    y_corner2 = y_line(idx_distances(max(1,end-corner_tolerance)));

                    % Make sure that the segment is longer than strel radius
                    % (i.e. no false positive)
                    if sqrt((x_corner1-x_corner2)^2+(y_corner1-y_corner2)^2) > (size(strel_overlap.Neighborhood,1)+1)/2
                        line_corners{i}(y_corner1,x_corner1) = 1;
                        line_corners{i}(y_corner2,x_corner2) = 1;
                    else
                        % revoke the segment's "good" status
                        good_segments_idx(j) = NaN;
                    end
                else
        %             do nothing
                end
            else
                % do nothing
            end
        end

        % CORNER REMOVAL
        % For blobs composed of "good segments" bridged by corners/intersections, 
        % only keep the lowest- and highest-numbered "corners" touching each blob.
        good_segments_idx = good_segments_idx(~isnan(good_segments_idx));
        good_segments     = ismember(segments,good_segments_idx);
        bridged_blobs = bwlabel(logical(good_segments + mask_intersections + line_corners{i}) .* line_map,8);
        n_bridged_blobs = max(bridged_blobs,[],'all');
        good_corner_blobs = bwlabel(line_corners{i},8);     % These should just be individual dots along the line
        for j = 1:n_bridged_blobs
            bridged_blob = bridged_blobs==j;
            bridged_blob_corners = bridged_blob .* good_corner_blobs;

            if max(bridged_blob_corners,[],'all') > 0
                % If there are corners on the bridged blob, then add the
                % furthest-away corners to the "final corners" array.
                corner_values_sorted = unique(bridged_blob_corners);    % 1D array of values in `bridged_blob_corners`.
                min_corner = corner_values_sorted(2);   % Excludes 0.
                max_corner = corner_values_sorted(end);

                corners_map = corners_map ...
                             + (bridged_blob_corners==min_corner) ...
                             + (bridged_blob_corners==max_corner);
            end
        end
    end

    corners_map = logical(corners_map);
%     figure; imagesc(edges + edges_s*0.4 + 5*corners_map); axis image; colorbar
end

function [x_output,y_output] = poly_sortcorners(corners_map,edges)
%POLY_SORTCORNERS() takes a 2D map indicating n-gon corners,
% as outputted by poly_findcorners(), and returns their
% coordinates in an n-by-2 array ordered by connectivity.
% (e.g. the first point is connected to the second point,
% which is connected to the third, etc., before looping
% back to the first.)
% 
% Returns:
% --------
% x_output : double (n-by-1)
%   x-coordinates for the n vertices.
% y-output : double (n-by-1)
%   y-coordinates for the n vertices.

    edge_fatten_radius = round(length(edges)/55);       % Just wide enough to include jaggies.
    strel_corner = strel('disk',5);                     % Enlarges corners; any touching will be merged.
    strel_edge   = strel('disk',edge_fatten_radius);    % Enlarges `edges`; detects whether corners are truly connected
    edges_fat = imdilate(edges,strel_edge);
    edges_veryfat = imdilate(edges,strel('disk',edge_fatten_radius*3));

    % Merge corners that are too close
    corners_map_fat = imdilate(corners_map,strel_corner);
    n_corners     = max(bwlabel(corners_map),[],'all');
    n_corners_fat = max(bwlabel(corners_map_fat),[],'all');
    if n_corners > n_corners_fat
        corners_map_backup = corners_map;
        corners_map = zeros(size(edges));
        for i = 1:n_corners_fat
            blob = bwlabel(corners_map_fat)==i;
            stats = regionprops(blob);
            centroid = stats.Centroid;
            x_cen = int16(centroid(1));
            y_cen = int16(centroid(2));
            corners_map(y_cen,x_cen) = 1;
        end
    end
    % Merge multiple-pixel "corners"
    corners_labeled = bwlabel(corners_map);
    n_corners       = max(corners_labeled,[],'all');
    for i = 1:n_corners
        blob = corners_labeled==i;
        if sum(blob,'all') > 1
            corners_map(blob) = 0;
            stats = regionprops(blob);
            centroid = stats.Centroid;
            x_cen = int16(centroid(1));
            y_cen = int16(centroid(2));
            corners_map(y_cen,x_cen) = 1; 
        end
    end
    % Remove corners that are too far away from edge
    corners_map = corners_map .* edges_veryfat;

    % Begin looping through corners. 
    % If algorithm fails to connect all of them (e.g. because
    % `overlap_threshold` is too strict), lower it and try again!
    overlap_threshold = 0.95;
    while true
        % Label the corners
        corners_labeled = bwlabel(corners_map);
        n_corners     = max(corners_labeled,[],'all');
        corner_candidates = 1:n_corners;                    % Tracks yet-unconnected corners
        corner_order      = NaN(size(corner_candidates));   % Will indicate connectivity order

        active_corner_bw = 1;   % Active corner of `corners_labeled`
        active_corner    = 1;   % Order of active corner in final n-gon (will go 1, 2, 3, etc)

        while true
            corner_candidates(corner_candidates==active_corner_bw) = [];    % Remove active corner from candidates
            corner_order(active_corner) = active_corner_bw;
            if isempty(corner_candidates)
                % If no corners are remaining, the active corner will automatically
                % connect with corner number 1
                break
            end
            % Find distances to candidate corners, and check 
            corner_distances = NaN(size(corner_candidates));    % Distance from corner to active corner.
            corner_isedge    = zeros(size(corner_candidates));  % Whether corner and active corner are connected by `edges`.
            corner_overlaps  = zeros(size(corner_candidates));  % Overlap fractions for each corner.
            idx = find(corners_labeled==active_corner_bw);
            [y_cornerA,x_cornerA] = ind2sub(size(edges),idx);      % x- and y-coordinates of edge.
            for i=corner_candidates
                idx = find(corners_labeled==i);
                [y_cornerB,x_cornerB] = ind2sub(size(edges),idx);      % x- and y-coordinates of edge.
                corner_distances(corner_candidates==i) = sqrt((x_cornerA-x_cornerB)^2+(y_cornerA-y_cornerB)^2);
                % Connected?
                [x_b,y_b] = bresenham(x_cornerA,y_cornerA,x_cornerB,y_cornerB);
                line_len = length(x_b);
                line_map = zeros(size(edges));          % Will indicate connecting line between corners
                for j = 1:line_len
                    x = x_b(j);
                    y = y_b(j);
                    line_map(y,x) = line_map(y,x)+1;
                end
        %             if active_corner_bw==6
        %                 figure; imagesc(line_map  + 2*corners_map + 0.2*edges_fat); axis image; colorbar;
        %                 width=1366;height=768;
        %                 set(gcf,'position',[200,200,width,height]);
        %             end
                overlap_fraction = abs(sum(line_map .* edges_fat,'all') / sum(line_map,'all'));
                corner_overlaps(corner_candidates==i) = overlap_fraction;
                if overlap_fraction > 0.95
                    % They're connected by an edge!
                    corner_isedge(corner_candidates==i) = 1;
                else
                    % They're not connected by an edge.
                    % Remove corner B from consideration.
                    corner_isedge(corner_candidates==i) = 0;
                end
            end

            % Workaround, if none of them turn out to be attached: pick the one
            %   that's most likely to be attached by edge!
            if ~ismember(1,corner_isedge)
                overlap_score = corner_overlaps ./ corner_distances;     % Higher overlap & shorter distances are better.
                max_score = max(overlap_score);
                corner_isedge(overlap_score==max_score) = 1;
%                 active_corner_bw
%                 corner_candidates
%                 corner_distances
%                 corner_overlaps
%                 overlap_score
%                 corner_order
%                 warning('poly_sortcorners() : Could not find corner connected by edge; least bad corner selected.');
            end

            % Find shortest distance to connected corner
            min_distance = min(corner_distances(logical(corner_isedge)));
            next_corner  = min(corner_candidates(corner_distances==min_distance));  % The "min()" is in case there are equidistant corners.

            % `active_corner` is attached to `next_corner`!
            active_corner_bw = next_corner;
            active_corner    = active_corner+1;
        end

        if sum(isnan(corner_order)) > 0
            overlap_threshold = overlap_threshold - 0.05;
            if overlap_threshold <= 0.60
                % If corners stubbornly fail to connect to each other:
                corner_candidates
                corner_distances
                corner_order
                error("poly_sortcorners() : Corners stubbornly failing to connect. Investigate.");
                % (!!!) Possible solution if this happens:
                % Either increase edges_fatten_radius to be more accomodating,
                % or decrease edges_fatten_radius so "bad" corners are left out
                % entirely.
            end
        else
            % If all corners are connected, exit the loop!
            break
        end
    end

    % Loop through corners in specified `corner_order`, and return coordinates!
    x_output = NaN(size(corner_order));
    y_output = NaN(size(corner_order));
    if length(corner_order)>1
        for ii=1:length(corner_order)
            idx = find(corners_labeled==corner_order(ii));
            [y_output(ii),x_output(ii)] = ind2sub(size(edges),idx);      % x- and y-coordinates of edge.
        end
    else
        % Return bogus "placeholder" corners that will result in a hilariously low score.
        disp("poly_sortcorners() : Corner sorting failed.");
        x_output = [1, 1, size(edges,2), size(edges,2)];
        y_output = [1, size(edges,1), size(edges,1), 1];
    end
end

function corners = poly_cleancorners(corners,edges)
%POLY_CLEANCORNERS() simply removes any corners from the final n-by-2 array
% that fall directly in between their adjacent corners (i.e. lie directly
% on a straight line).    
    % Dilate corners to see if they're intersected by lines connecting the
    % previous and next corners. If so, it will be removed.
    corner_radius = 3;
    bad_corners = [];   % List of corners that will be removed.
    
    % Loop through corners; check the one in the middle.
    n_corners = size(corners,1);
    for j = 1:n_corners
        % corner_current
        x_current = corners(j,1);
        y_current = corners(j,2);
        % corner_middle
        x_mid = corners(rem(j+0,n_corners)+1,1);
        y_mid = corners(rem(j+0,n_corners)+1,2);
        % corner_end
        x_end = corners(rem(j+1,n_corners)+1,1);
        y_end = corners(rem(j+1,n_corners)+1,2);
        
        % Draw the (fattened) middle corner on the map.
        corner_mid_map = zeros(size(edges));
        corner_mid_map(y_mid,x_mid) = 1;
        corner_mid_map = imdilate(corner_mid_map,strel('disk',corner_radius));
        
        % Draw a line from corner_current to corner_end.
        % If corner_mid (fattened) is on that line, flag it as bad!
        [x_b,y_b] = bresenham(x_current,y_current,x_end,y_end);
        line_len = length(x_b);
        for k = 1:line_len
            x = x_b(k);
            y = y_b(k);
            corner_mid_map(y,x) = corner_mid_map(y,x) + 1;
        end
        if max(corner_mid_map,[],'all') > 1
            % Middle corner needs to be removed!
            bad_corners = [bad_corners,rem(j+0,n_corners)+1];
        end        
    end
    corners(bad_corners,:) = [];
end

function corners = poly_untangle(corners,edges,loop)
% Identifies where lines drawn between adjacent polygon corners cross each
% other, and attempts to "untangle" the lines accordingly.

    % Initialize
    strel_widen = strel('disk',3);
    corners_backup = corners;

    % Draw polygon traced by vertices, one line at a time
    n_corners = size(corners,1);
    corner_b  = cell(n_corners,1);      % x- and y-coords for each (j->jp1) line
    map_polygon = zeros(size(edges));
    for j = 1:n_corners
        jp1 = rem(j+0,n_corners)+1;  % Next corner
        x_current = corners(j,1);
        y_current = corners(j,2);
        % corner_next
        x_next = corners(jp1,1);
        y_next = corners(jp1,2);
        
        % Draw a line from current corner to next corner!
        [x_b,y_b] = bresenham(x_current,y_current,x_next,y_next);
        x_b(1) = [];    y_b(1) = [];    % Ignore first point.
        corner_b{j} = [x_b,y_b];        % Save Bresenham line for each (j->jp1) line
        line_len = length(x_b);
        map_polygon_backup = map_polygon;   % Backup, without the newest line drawn.
        for k = 1:line_len
            x = x_b(k);
            y = y_b(k);
            map_polygon(y,x) = map_polygon(y,x) + 1;
        end
        
        % Check for intersection
        intersection_found = max(map_polygon,[],'all') > 1;
        map_intersection = map_polygon > 1;     % 2D; "1" wherever intersection exists
        pix_skip = 5;   % Number of pixels being ignored when using "fat" line
        if ~intersection_found && line_len>(pix_skip*2)
            % Check for "fat intersection". Requires line length of 10
            % or more (since that's how many pixels will be removed,
            % just in case).
            x_b_f = x_b;
            y_b_f = y_b;
            x_b_f(1:pix_skip) = [];    x_b_f(end-pix_skip+1:end) = [];
            y_b_f(1:pix_skip) = [];    y_b_f(end-pix_skip+1:end) = [];
            line_len = length(x_b_f);
            map_shortline = zeros(size(map_polygon));    % 2D map of line, >`pix_skip` pixels away from start & end.
            for k = 1:line_len
                x = x_b_f(k);
                y = y_b_f(k);
                map_shortline(y,x) = map_shortline(y,x) + 1;
            end
            map_shortline_fat = imdilate(map_shortline,strel('disk',1));
            map_intersection = map_shortline_fat .* map_polygon_backup;
            intersection_found = max(map_intersection,[],'all') > 0;
            if intersection_found
%                 disp('FAT INTERSECTION FOUND!');
            end
        end
        
        % Untangle, if the intersection is serious enough!
        if intersection_found
            idx = find(map_intersection);
            idx = idx(1);   % Only consider 1 intersection at a time.
            [y_int,x_int] = ind2sub(size(edges),idx);
            
            % Find which existing line is being intersected
            for k = 1:n_corners
                if ismember([x_int,y_int], corner_b{k}, 'rows')
                    break
                end
            end
            kp1 = rem(k+0,n_corners)+1;
%             disp("TANGLE FOUND between lines "+j+"-"+jp1+" and "+k+"-"+kp1);
            
            % A tangle is found between j-jp1 and k-kp1.
            % Typically, this means that all corners between (inclusive)
            %  j-kp1 OR jp1-k (whichever pair is closer) have their order
            %  reversed.
            % EXCEPTION: If j==kp1 or jp1==k, then we'll need to find out
            %  which of the two lines is experiencing more overlap.
            %  e.g. Consider two lines (corners labeled) that look like:
            %       4----6==========5
            %  ^ Here, j=5, jp1=6,
            %          k=4, kp1=5.
            %  Since line 5-6 is more overlapped than 4-5, we would switch
            %  around corners 5 and 6.
            if (j==kp1 || jp1==k) && (sum(map_intersection,'all') > 2)
                % Identify which of the corners appears twice.
                joint = [j jp1 k kp1];  % Will be the corner that appears twice.
                [other,idx] = unique(joint);
                joint(idx) = [];
                other(other==joint) = [];
                overlap_score = zeros(size(other));
                % Loop through the "other" corners; see which lines are
                % mostly overlapped!
                for jj=other
                    line_map = zeros(size(edges));
                    corner = min(jj,joint);
                    x_b = corner_b{corner}(:,1);
                    y_b = corner_b{corner}(:,2);
                    
                    % Draw a line!
                    line_len = length(x_b);
                    for k = 1:line_len
                        x = x_b(k);
                        y = y_b(k);
                        line_map(y,x) = 1;
                    end
                    line_wide = imdilate(line_map,strel_widen);
                    % Calculate overlap! (Higher score means more overlap.)
                    overlap_score(other==jj) = sum(line_wide .* map_polygon,'all') / sum(line_map,'all');
                end
                % Whichever "other" corner has more overlap gets swapped
                % with "joint"!
                swap1 = other(overlap_score==max(overlap_score));
                swap2 = joint;
                corners([swap1,swap2],:) = corners([swap2,swap1],:);
%                 disp("SWAPPED "+swap1+" and "+swap2)
                break
            elseif j~=kp1 && kp1~=k
                % NOTES:
                % - j>k
                % - kp1 is 1 closer to j than k is
                % We want to find distance between j-kp1, and jp1-k.
                % Whichever is closer, everything between that pair
                % (increasing from k~ to j~)
                dist_jkp1 = corner_distance(j,kp1,n_corners);
                dist_jp1k = corner_distance(jp1,k,n_corners);
                if dist_jkp1 < dist_jp1k
                    corners_to_reverse = corners_listbetween(kp1,j,n_corners);
                else
                    corners_to_reverse = corners_listbetween(jp1,k,n_corners);
                end
                % Reverse the corners!
%                 disp("REVERSED the order of ["+num2str(corners_to_reverse)+"]")
                corners(corners_to_reverse,:) = corners(flip(corners_to_reverse),:);
                break
            else
                % If the intersection isn't enough to be worth swapping for
                map_polygon(map_polygon > 1) = 1;
%                 disp("Must be the wind.");
                continue
            end
        end
    end
    if ~isequal(corners,corners_backup) && (loop<=n_corners)
        % Try again and see if more can be untangled!
%         disp('LOOP');
        corners = poly_untangle(corners,edges,loop+1);
    elseif loop>=n_corners
        % It's no good. Return bogus corners that will result in a hilariously low score.
        disp("poly_untangle() : Edge 'untangling' ran into an infinite loop. Canceling.");
%         x_placeholder = [1, 1, size(edges,2), size(edges,2)];
%         y_placeholder = [1, size(edges,1), size(edges,1), 1];
%         corners = [x_placeholder',y_placeholder'];
        corners = corners;
    else
        % Do nothing; let it end here.
    end
end

function itbuffer = createLineIterator(P1, P2, im)
%     Parameters:
%         -P1: a numpy array that consists of the coordinate of the first point (x,y)
%         -P2: a numpy array that consists of the coordinate of the second point (x,y)
%         -im: the image being processed
%     Returns:
%         -it: a numpy array that consists of the coordinates and intensities of each pixel in the radii (shape: [numPixels, 3], row = [x,y,intensity])
    imageH = size(im,1);
    imageW = size(im,2);
    P1X = double(P1(1));
    P1Y = double(P1(2));
    P2X = double(P2(1));
    P2Y = double(P2(2));
    % difference and absolute difference between points
    % used to calculate slope and relative location between points
    dX = P2X - P1X;
    dY = P2Y - P1Y;
    dXa = abs(dX);
    dYa = abs(dY);
    % predefine numpy array for output based on distance between points
    %     itbuffer = np.empty(shape=(np.maximum(dYa, dXa), 3), dtype=np.float32)
    %     itbuffer.fill(np.nan)
    itbuffer = nan(max(dYa, dXa), 3, 'double');

    % Obtain coordinates along the line using a form of Bresenham's algorithm
    negY = P1Y > P2Y;
    negX = P1X > P2X;
    if P1X == P2X  % vertical line segment
        itbuffer(:,1) = P1X;
        if negY
    %             itbuffer(:, 2) = np.arange(P1Y - 1, P1Y - dYa - 1, -1);
            itbuffer(:, 2) = [(P1Y - 1):-1:(P1Y - dYa)];
        else
    %             itbuffer(:, 2) = np.arange(P1Y + 1, P1Y + dYa + 1);
            itbuffer(:, 2) = [(P1Y + 1):1:(P1Y + dYa)];
        end
    elseif P1Y == P2Y  % horizontal line segment
        itbuffer(:, 2) = P1Y;
        if negX
    %             itbuffer(:, 1) = np.arange(P1X - 1, P1X - dXa - 1, -1);
            itbuffer(:, 1) = [(P1X - 1):-1:(P1X - dXa)];
        else
    %             itbuffer(:, 1) = np.arange(P1X + 1, P1X + dXa + 1);
            itbuffer(:, 1) = [(P1X + 1):1:(P1X + dXa)];
        end
    else  % diagonal line segment
        steepSlope = dYa > dXa;
        if steepSlope
            slope = double(dX) / double(dY);
            if negY
    %                 itbuffer(:, 2) = np.arange(P1Y - 1, P1Y - dYa - 1, -1)
                itbuffer(:, 2) = [(P1Y - 1):-1:(P1Y - dYa)];
            else
    %                 itbuffer(:, 2) = np.arange(P1Y + 1, P1Y + dYa + 1);
                itbuffer(:, 2) = [(P1Y + 1):1:(P1Y + dYa)];
            end
            itbuffer(:, 1) = [int16(slope * (itbuffer(:, 2) - double(P1Y)) + double(P1X))];
        else
            slope = double(dY) / double(dX);
            if negX
    %                 itbuffer(:, 1) = np.arange(P1X - 1, P1X - dXa - 1, -1)
                itbuffer(:, 1) = [(P1X - 1):-1:(P1X - dXa)];
            else
    %                 itbuffer(:, 1) = np.arange(P1X + 1, P1X + dXa + 1)
                itbuffer(:, 1) = [(P1X + 1):1:(P1X + dXa)];
            end
            itbuffer(:, 2) = [int16(slope * (itbuffer(:, 1) - double(P1X)) + double(P1Y))];
        end
    end
    % Remove points outside of image
    colX = itbuffer(:, 1);
    colY = itbuffer(:, 2);
    idx_good = (colX>=1) & (colY >=1) & (colX <= imageW) & (colY <= imageH);
    colX = colX(idx_good);
    colY = colY(idx_good);
    itbuffer = [colX,colY,nan(size(colX,1),1)];
    % itbuffer = itbuffer((colX >= 0) & (colY >= 0) & (colX < imageW) & (colY < imageH));
    % % Get intensities from img ndarray
    % itbuffer(:, 3) = im(uint16(itbuffer(:, 2)), uint16(itbuffer(:, 1)));
    idx_buffer = sub2ind(size(im),itbuffer(:, 2), itbuffer(:, 1));
    itbuffer(:, 3) = im(idx_buffer);
end