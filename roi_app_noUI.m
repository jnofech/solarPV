function roi_app_noUI(name, hmap, map, pix2m, isroof, outputname, ~)
% Basically `roi_app_test`, but without a UI that pops up.
% This is a workaround to resolve a bug in MATLAB's App Designer where an
% app's UI somehow opens _after_ running and sending a close request.
            
    % Get map resolution
    width = size(hmap,2);
    height = size(hmap,1);
    ScreenSize = get(groot,'ScreenSize');
    screen_height = ScreenSize(4);
    screen_width  = ScreenSize(3);

    % Define variables
    hmap = hmap - min(hmap,[],'all');   % Lowest value must be zero.

    % Generate "base mask"
    % Dilate `isroof`.
    se = strel('disk',ceil(2 / pix2m));                 % Rooftops in "base mask" must be no more than 2m apart
    isroof_fat = imdilate(isroof,se);
    isroof_fat = imfill(isroof_fat,'holes');

    % Draw a circle, and "select" any surface blobs in contact with it.
    % If circle isn't touching anything, enlarge it until it is.
    map_dot = zeros(size(hmap));
    map_dot(round(end/2),round(end/2)) = 1;     % Put a tiny dot in the middle.
    circ_radius = ceil(min(size(hmap))/10);     % This dot will turn into a circle of r=`circ_radius`.
    while true
        circ_radius = circ_radius + 20;
        se = strel('disk',circ_radius);
        map_circle = imdilate(map_dot,se);
        if max(map_circle + isroof_fat,[],'all') > 1
            break
        end
    end

    % Create "base mask" (i.e. the "default" mask without modifications)!
    blobs = bwlabel(isroof_fat,4);
    n_blobs = max(blobs,[],'all');
    base_mask = zeros(size(hmap));

    for i = 1:n_blobs
        blob = (blobs==i);
        % Check if blob is intersecting circle.
        % If so, add that blob to base mask!
        if max(blob + map_circle,[],'all') > 1
            base_mask = base_mask + blob;
        end
    end
    % Initialize hmap display with this `base_mask`
    hmap_disp = hmap;
    hmap_disp = hmap_disp + ~base_mask*max(hmap,[],'all');  % "Mask" from base_mask.
    hmap_disp(~isfinite(hmap_disp)) = max(hmap,[],'all');   % "Mask" out any NaN regions.
    hmap_disp = hmap_disp/(2*max(hmap,[],'all'));           % "Normalize", such that final display uses the combined cmap.

    % Save file and exit
    mask = base_mask;
    mask_id = randi([1,999999999],1);   % Randomly-selected "mask ID". This is so that polygon approximation is performed per mask.
    if endsWith(outputname,".mat")
        app.outputname = outputname;
    else
        app.outputname = outputname + ".mat";
    end
    save(app.outputname, 'mask', 'mask_id');
    disp("Saved output mask to "+app.outputname+". Closing.");
end