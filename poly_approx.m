function table_polygons = poly_approx(name,info,roofs,roofs_p,mask_id,pix2m,kmeans_mode,force_run,output_path_save)
%POLY_APPROX generates approximate polygons for all of a building's roofs,
%  returning them (and their relative confidence) as a table. Note that if
%  polygons already exist for a corresponding building mask, those results
%  will be reused.
    % Output 
    outputname = output_path_save+name+"_polygons_"+mask_id+".mat";

    if ~isfile(outputname) || force_run==true
        disp(name+": Beginning polygon approximation for mask #"+mask_id+".");
        n_rows = length(info.roof);
        polygons            = cell(n_rows,1);   % Corners [x,y], in order of connectivity.
        polygons_confidence = zeros(n_rows,1);  % Confidence (from 0 to 1-ish) of said corners.
        polygons_obstacles  = cell(n_rows,1);   % Corners [x,y] of OBSTACLES, in order of connectivity.
        polygons_combined   = cell(n_rows,1);   % Corners [x,y] of ROOFTOP+OBSTACLES (i.e. the same information in `polygons` and `polygons_obstacles`),
                                                % separated by NaNs.
        for i = 1:n_rows
            disp("Roof "+i+" of "+n_rows)
            if info.type(i) ~= "irregular" && (~info.toosmall(i))
                [rows, cols] = find(roofs_p==i);
                row1 = min(rows);
                col1 = min(cols);
                [corners,confidence] = find_poly(roofs==i,roofs_p==i,kmeans_mode);  % Returns corners in "whole building" pixels. Clockwise.
                corners_obstacles    = find_poly_obstacles(roofs==i,corners,pix2m); % Returns obstacle corners in "whole building" pixels. Counterclockwise.
                corners(:,1) = corners(:,1) - col1;
                corners(:,2) = corners(:,2) - row1;
                for ii = 1:size(corners_obstacles,1)
                    corners_obstacles{ii}(:,1) = corners_obstacles{ii}(:,1) - col1;
                    corners_obstacles{ii}(:,2) = corners_obstacles{ii}(:,2) - row1;
                end
                % ^ Same coordinates as the binary itself! Note that
                % some negative values will appear.
                polygons{i}             = corners;
                polygons_confidence(i)  = confidence;
                polygons_obstacles{i}   = corners_obstacles;
                if isempty(polygons_obstacles{i})
                    polygons_obstacles{i}   = [NaN, NaN];
                end
                
                % Combine rooftop poly's with obstacle poly's, as per Nima's request
                polygons_combined{i} = polygons{i};
                if isequaln(polygons_obstacles{i},[NaN,NaN])
                    % do nothing
                else
                    obstacles = polygons_obstacles{i};
                    for ii = 1:size(polygons_obstacles{i},1)
                        obstacle = obstacles{ii};
                        polygons_combined{i} = [polygons_combined{i}; [NaN,NaN]; obstacle];
                    end
                end

            else
                polygons{i}             = [NaN, NaN];
                polygons_confidence(i)  = 0;
                polygons_obstacles{i}   = [NaN, NaN];
            end
        end
        table_polygons = table(polygons, polygons_confidence, polygons_obstacles, polygons_combined);
        save(outputname,'table_polygons')
    else
        disp(name+": Reading polygons from mask ID#"+mask_id+".");
    end
    load(outputname,'table_polygons');
%     disp("Poly table height: "+height(table_polygons));
end

