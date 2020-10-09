function map = map_read(name,input_path_photos,output_path_save)
%MAP_READ returns a 2D satellite photo for a building `name`, in pixels.
%Will read '.png' and '.jpg' files.
%
%   Parameters:
%   -----------
%   name : string
%       Name of building. Case-sensitive.
%
%   Outputs:
%   --------
%   `map = map_read(name)`:
%       `map` is a colour image (3D array) of the building.
    name = string(name);
%     input_path_photos = "images/";
%     output_path_save = "output/";
    fname  = input_path_photos+name;
    output_name = output_path_save+name+"_map_rotation.mat";
    
    % Read file
    if isfile(fname+".png")
        map = imread(fname+".png");
    elseif isfile(fname+".jpg")
        map = imread(fname+".jpg");
    else
        map = [];
        warning(fname+".png and "+fname+".jpg do not exist.");
        return;
    end
    
    if ~isfile(output_name)
        % Get direction of map.
        figure; 
        width=1366;height=768;
        set(gcf,'position',[200,200,width,height]);imagesc(map); axis image; 
        title(fname);
        cmap_dir = 'None';
        while strcmp(cmap_dir,'None')
            button = questdlg('Select '+name+' colourmap direction. (Press ESC to skip this step)', ...
                              'FUNCTION(): '+name+' colourmap direction', 'North', 'East', 'Next', 'North');
            if strcmp(button, 'North') | strcmp(button, 'East') | strcmp(button,'')
                cmap_dir = button;
            else
                button = questdlg('Select '+name+' colourmap direction. (Press ESC to skip this step)', ...
                                  'FUNCTION(): '+name+' colourmap direction', 'South', 'West', 'Previous', 'South');
                if strcmp(button, 'South') | strcmp(button, 'West') | strcmp(button,'')
                    cmap_dir = button;
                elseif strcmp(button, 'Previous')
                    cmap_dir = 'None';
                end
            end
        end
        close;
        if strcmp(button,'')
            warning('No colourmap direction selected. Assuming North.');
            button = 'North';
        end
        % Rotate the map North.
        switch button
            case 'North'
                % Do nothing
                theta_CCW = 0;
            case 'East'
                % Rotate 90deg CW
                theta_CCW = 270;
            case 'South'
                % Rotate 180deg CW
                theta_CCW = 180;
            case 'West'
                % Rotate 90deg CCW
                theta_CCW = 90;
            otherwise
                error("Invalid `button` output")
        end
        % Save rotated angle
        save(output_name,'theta_CCW');
        % Rotate map
        map = imrotate(map,theta_CCW,'bilinear');   % Rotates map, filling gaps with zeros.
%         imagesc(map); axis image;
%         title("Rotated colour map. If incorrect, delete '"+replace(output_name,'_','\_')+"' and retry!")
    else
        % Load rotated angle
        load(output_name,'theta_CCW');
        % Rotate map
        map = imrotate(map,theta_CCW,'bilinear');   % Rotates map, filling gaps with zeros.
    end
    
end