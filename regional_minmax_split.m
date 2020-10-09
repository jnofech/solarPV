function [isroof_modified] = regional_minmax_split(Z,dZdX,dZdY,roofs,roof_isbad,pix2m)
% REGIONAL_MINMAX_SPLIT() detects maximum min/max values along a line of
% ~4m length, for angles in steps of 45 degrees. Regional min/max locations
% are identified as "ridges" and can split blobs into smaller sub-blobs
% based on differences in slope and angle.
% This code essentially prepares _gabled/butterfly rooftops_ for use in the
% Building() object.
    isroof_modified = zeros(size(roofs));
    n_blobs = max(roofs,[],'all');
    for i=1:n_blobs
        slice_full = (roofs==i);                % Binary array of each blob
        [rows, cols] = find(slice_full);
        row1 = min(rows); row2 = max(rows);
        col1 = min(cols); col2 = max(cols);

        slice = slice_full(row1:row2,col1:col2);        % Cropped binary array containing each blob
        slice_mask = double(slice);
        slice_mask(slice==0) = nan;                     % nan where roof isn't
        slice_isgood = (~roof_isbad(row1:row2,col1:col2)).*slice_mask;
        slice_Z = Z(row1:row2,col1:col2).*slice_mask;   % Heights of blob
        slice_Zx = dZdX(row1:row2,col1:col2).*slice_mask;   % Slopes of blob
        slice_Zy = dZdY(row1:row2,col1:col2).*slice_mask;

        % Find regional min/max w/ linear structuring elements (i.e. shows 
        % minimum value within a line centered at each pixel).
        % For each of min/max, use a LINE-SHAPED strel of ~2.5m radius, 
        % oriented in X, Y, and diagonal directions.
        %   If a point is the regional maximum or minimum within ANY of these
        %   directions, then it may be a ridge pointing upwards or downwards!
        se_length = ceil(5/pix2m);
        slice_modified = slice;                 % Parts will be cut out of this
        slice_mask     = slice*0;               % (will be) final mask applied to slice

        for mult = [1,-1]           % For [min,max] regions
            map = mult*slice_Z;
            map(isnan(map)) = 99999;        % Remove NaNs from consideration.
            for se_angle = [0,90,45,135]    % -, |, /, \ directions
    %         for se_angle = [0,90]           % -, | directions
    %         for se_angle = 90
                se = strel('line',se_length,se_angle);
                masks = (slice_Z == mult*imerode(map,se));  % Erode for mins
                masks = imdilate(masks,strel('disk',1));    % Thicken mask by 1
                    % EXPERIMENTAL FIX FOR REAL RIDGES NOT SPLITTING SLICE
                    masks = imdilate(masks,strel('disk',1));
                    masks = imerode(masks,strel('disk',1));
                % Consider ONLY the parts of mask that split the slice into 
                % multiple distinct blobs!
                mask_miniblobs = bwlabel(masks,4);
    %             n_mask_miniblobs = max(mask_miniblobs,[],'all');
    %             for k = 1:n_mask_miniblobs
    %                 mask_miniblob = mask_miniblobs==k;
                    % NOTE:
                    % It's "safer" to loop through one miniblob at a time, as
                    % mask blobs that overlap may cause issues. However, this
                    % is rare, as the process is done for one `se_angle` at a
                    % time, so they're unlikely to cross in the first place.
                    mask_miniblob = masks;
    %                 figure; imagesc(slice + 5*mask_miniblob); axis image; colorbar;
                [slice_split,mask_split] = blob_split_matchslope(slice,slice_isgood,slice_Zx,slice_Zy,mask_miniblob);
    %                 figure; imagesc(slice_split + 3*mask_split); axis image; colorbar;
                    slice_modified = slice_modified & slice_split;
                    slice_mask = slice_mask + mask_split;
    %             end
            end

            % (!!!!!) MAYBE... Split the slice, and then LOOP THROUGH ALL OF
            % THE `regional` BLOBS AGAIN and see if another split happens! Or
            % just, like... finish getting the split blobs, and then run the
            % whole code one more time.
            % If this works, this may provide a solution for hip roofs!
        end

        % REPAIR ROOFTOPS! For the above, miscalculations may occur because of
        % miniblobs that SHOULD be separate, but AREN'T due to being separated
        % by a differently-angled strel (or peak/trough).
        slice_modified = bwareaopen(slice_modified, min([round(nansum(slice,'all')*0.025),40]), 4); % Remove small chunks.
        slice_mask = slice_mask & slice;
            % EXPERIMENTAL FIX FOR MULTIPLE CUTS "MEETING" AT A CORNER
            % Erode (slice_modified + slice_mask) by 3 pixels.
            % Get the overlap of modified_slices (dilated by 5pix or somethin).
            % Wherever the overlap is 3 or more, multiply the mask by the
            % eroded thingy above. (Eats away at wherever mask is touching
            % background IF that region has 3 or more surfaces nearby.)
            mask_mask = imerode(slice_modified+slice_mask,strel('disk',4));
            % Dilate the split slice by a few pixels. Wherever the thickened
            % regions overlap 3 or more times is where the mask might have to
            % be cut away a bit.
            overlap = zeros(size(slice));      % (will be) replacement for `mask`
            modified_slices = bwlabel(slice_modified,4);
            for i = 1:max(modified_slices,[],'all')
                miniblob = (modified_slices==i);
                overlap = overlap + (imdilate(miniblob,strel('disk',6)) - miniblob);
                % ^ Dilation reaches out <~3 pixels from miniblob
            end
            slice_mask = slice_mask .* ~((overlap >= 3) & ~mask_mask);

        slice_repaired = blob_repair(slice_modified,slice_mask,slice_isgood,slice_Zx,slice_Zy);

%         isroof_modified(row1:row2,col1:col2) = isroof_modified(row1:row2,col1:col2) | slice_modified;
        isroof_modified(row1:row2,col1:col2) = isroof_modified(row1:row2,col1:col2) | slice_repaired;
    end 
end

function [slice_repaired] = blob_repair(slice_split,mask,slice_isgood,slice_Zx,slice_Zy)
% Given a slice that had a mask applied, which may be separated into
% component blobs, "repair" any of the cuts where the bridged slices would
% behave identically afterwards.
    slope_minslant = 2*0.3/12;   

    slice_split_miniblobs = bwlabel(slice_split,4);
    n_split_miniblobs = max(slice_split_miniblobs,[],'all');
    bridges_final = zeros(size(slice_split));   % (will be) final "bridges" added to slice_split
    
    slice_Zx = slice_Zx .* slice_isgood;
    slice_Zy = slice_Zy .* slice_isgood;
    
    slice_mask_miniblobs = bwlabel(mask,4);
    n_mask_miniblobs = max(slice_mask_miniblobs,[],'all');
    % Loop through "potential bridges" of mask. Check if each is bridge
    for ii = 1:n_mask_miniblobs
        potential_bridge = (slice_mask_miniblobs == ii);
        miniblobs_potentially_bridged = bwlabel((slice_split + potential_bridge),4);
        n_miniblobs_potentially_bridged = max(miniblobs_potentially_bridged,[],'all');
        % If this potential bridge is actually a bridge...
        if n_miniblobs_potentially_bridged < n_split_miniblobs
            % The following code works best if the bridge only connects
            % 2 miniblobs, rather than 3 or more.
            if n_split_miniblobs - n_miniblobs_potentially_bridged >= 2
%                 error("What do");
            end
            slice_connectedbybridge = slice_split;      % (will be) slices connected by current "potential bridge" (labeled).
            % Ignore miniblobs that aren't linked by this bridge
            for i = 1:n_miniblobs_potentially_bridged
                miniblob = (miniblobs_potentially_bridged==i);
                % If miniblob isn't touching the "potential bridge", then
                % remove it from consideration.
                if max(bwlabel(miniblob+potential_bridge),[],'all')>1
                    slice_connectedbybridge = slice_connectedbybridge - miniblob;
                end
            end
            % Finalize miniblobs connected by bridge
            slice_connectedbybridge = bwlabel(logical(slice_connectedbybridge),4);
            n_slice_connectedbybridge = max(slice_connectedbybridge,[],'all');
            % Get information about each of these miniblobs
            miniblob_Zx_means = zeros(n_slice_connectedbybridge,1);
            miniblob_Zy_means = zeros(n_slice_connectedbybridge,1);
            miniblob_m_means  = zeros(n_slice_connectedbybridge,1);
            miniblob_mean_m   = zeros(n_slice_connectedbybridge,1);
            miniblob_directions = zeros(n_slice_connectedbybridge,1);
            miniblob_substantial= zeros(n_slice_connectedbybridge,1);
            miniblob_slanted  = zeros(n_slice_connectedbybridge,1);
            miniblob_rounded  = zeros(n_slice_connectedbybridge,1);
            for jj = 1:n_slice_connectedbybridge    % For each connected miniblob
                miniblob = (slice_connectedbybridge==jj);
                where_miniblob = find(miniblob);
                miniblob_Zx_means(jj) = nanmean(slice_Zx(where_miniblob),'all');
                miniblob_Zy_means(jj) = nanmean(slice_Zy(where_miniblob),'all');
                miniblob_mean_m(jj)  = nanmean(sqrt(slice_Zx(where_miniblob).^2 + slice_Zy(where_miniblob).^2));     % Means of magnitude of slope. (e.g. Sphere treated as "slanted" with direction implied to be reliable.)
                miniblob_m_means(jj)  = sqrt(nanmean(slice_Zx(where_miniblob),'all').^2 + nanmean(slice_Zy(where_miniblob),'all').^2);  % Magnitude of means of slope. (e.g. Sphere treated as "flat", so direction is unreliable.)
                miniblob_directions(jj) = meanangle(-angle(slice_Zy(where_miniblob) + 1i*slice_Zx(where_miniblob)) * 180/pi,'all'); % Mean of angles of blob, CCW from North.
%                 miniblob_directions(jj) = -angle(miniblob_Zy_means(jj) + 1i*miniblob_Zx_means(jj)) * 180/pi;            % Angle from mean slopes of blob, CCW from North.
                miniblob_substantial(jj)= sum(miniblob,'all') > (nansum(logical(slice_connectedbybridge),'all')*0.05); % If miniblob is >5% of total bridged region

                if miniblob_substantial(jj)
                    miniblob_slanted(jj)  = miniblob_m_means(jj) > slope_minslant;
                    miniblob_rounded(jj)  = miniblob_mean_m(jj) > (2.5*miniblob_m_means(jj)+0.05);
                    % (!!!) Still needs calibration!
                else
                    miniblob_slanted(jj)  = miniblob_m_means(jj) > 5*slope_minslant;
                    miniblob_rounded(jj)  = miniblob_mean_m(jj) > (3*miniblob_m_means(jj)+0.05);
                    % Small chunks have a less reliable slope, so they
                    % need to be steeper to be distinguished from the 
                    % rest of the slice.
                end
            end
            anglediffs = @(angles) abs(angles - angles') - 360*(abs(angles - angles')>180);
            % FINAL CRITERIA:
            if ((max(abs(anglediffs(miniblob_directions)),[],'all') < 25 ...
             && ~ismember(false,miniblob_slanted)) ...  
             || ~ismember(true,miniblob_slanted)) ...  % If miniblobs have similar slope angle OR are all flat...
             && (1.3*min(miniblob_m_means)+0.1) > max(miniblob_m_means) ... %&& ~ismember(false,miniblob_substantial) ...
             && (~ismember(true,miniblob_rounded)) % ... and slope mags are similar and none of the miniblobs are "rounded"
                bridges_final = bridges_final + potential_bridge;
            else
            end
        end
    end
    slice_repaired = slice_split + bridges_final;
end

function [slice_split,mask_split] = blob_split_matchslope(slice,slice_isgood,slice_Zx,slice_Zy,mask)
% Given a slice, its slopes, and a mask, this code will see if the mask
% separates the slice into multiple distinct regions. If so, it will try to
% split the slice into regions of distinct slope (magnitude and/or
% direction).
    slope_minslant = 0.3/12;   
    % Miniblobs must have this slope magnitude to be classified as
    % "slanted". If the blob is so small that their slope isn't reliable,
    % then they must be 5x steeper than this to be treated as "slanted".
    % For comparison, gentle rainfall slopes may be somewhere around ~0.05, 
    % or 0.6/12 (e.g. slanted roof of UAlberta's CCIS).
    
    slice_Zx = slice_Zx .* slice_isgood;
    slice_Zy = slice_Zy .* slice_isgood;
    % When calculating roof slopes, only consider `~_isbad` regions
    
    slice_split = slice .* ~mask;
    slice_split = bwareaopen(slice_split, round(nansum(slice,'all')*0.025), 4); % Remove small chunks.
    slice_split_miniblobs = bwlabel(slice_split,4);
    n_split_miniblobs = max(slice_split_miniblobs,[],'all');
    % Mask may involve "squiggly lines" that reach inside the slice without
    % splitting it. To remove these:
    % Dilate the remaining miniblobs by 3 pixels. Wherever the thickened
    % regions overlap is where the mask is placed.
    overlap = zeros(size(slice));      % (will be) replacement for `mask`
    for i = 1:n_split_miniblobs
        miniblob = (slice_split_miniblobs==i);
        overlap = overlap + (imdilate(miniblob,strel('disk',4)) - miniblob);
        % ^ Dilation reaches out <~3 pixels from miniblob
    end
        % EXPERIMENTAL FIX: 
        % To fix "gaps" appearing where mask was thicker than 3 pixels
        % (i.e. where miniblobs were too far appart for the 3-pixel overlap
        % to connect them):
        % `slice_split` is dilated by <~2pix, and the holes within
        % this "smoothed" slice identify the aforementioned "gaps".
        % These gaps are dilated by 1pix and added to the mask.
        slice_split_smooth = imdilate(slice_split,strel('disk',3));
        mask_gap = imdilate(imfill(slice_split_smooth,4,'holes') - slice_split_smooth,strel('disk',1));
    mask = logical((overlap >= 2) + mask_gap);
    mask = mask & slice;
    mask_split = zeros(size(mask));
    
%                 figure; imagesc(overlap); axis image; colorbar; title('Overlap');
%                 figure; imagesc(mask_gap); axis image; colorbar; title('Regions missed by overlap');
%                 figure; imagesc(slice + 3*mask_gap + 3*mask); axis image; colorbar; title('Improved mask (>=3)');
%                 figure; imagesc(slice + 5*mask); axis image; colorbar; title('Base slice');
%                 figure; imagesc(slice_split + 5*mask); axis image; colorbar; title('Final main mask');

    slice_split = slice .* ~mask;    % Final masked binary of slice
    slice_split_miniblobs = bwlabel(slice_split,4);
    n_split_miniblobs = max(slice_split_miniblobs,[],'all');
    
    % If a split has occurred, then only apply mask at the "bridge"!
    if n_split_miniblobs >= 2
        slice_bridge = zeros(size(slice_split)); % (will be) "bridge" part of mask
        slice_mask_miniblobs = bwlabel(mask,4);
        n_mask_miniblobs = max(slice_mask_miniblobs,[],'all');
        % Loop through "potential bridges" of mask. Check if each is bridge
        for ii = 1:n_mask_miniblobs
            potential_bridge = (slice_mask_miniblobs == ii);
            miniblobs_potentially_bridged = bwlabel((slice_split + potential_bridge),4);
            n_miniblobs_potentially_bridged = max(miniblobs_potentially_bridged,[],'all');
            % If this potential bridge is actually a bridge...
            if n_miniblobs_potentially_bridged < n_split_miniblobs
                % The following code works best if the bridge only connects
                % 2 miniblobs, rather than 3 or more.
                if n_split_miniblobs - n_miniblobs_potentially_bridged >= 2
%                     error("What do");
                end
                slice_connectedbybridge = slice_split;      % (will be) slices connected by current "potential bridge" (labeled).
                % Ignore miniblobs that aren't linked by this bridge
                for i = 1:n_miniblobs_potentially_bridged
                    miniblob = (miniblobs_potentially_bridged==i);
                    % If miniblob isn't touching the "potential bridge", then
                    % remove it from consideration.
                    if max(bwlabel(miniblob+potential_bridge),[],'all')>1
                        slice_connectedbybridge = slice_connectedbybridge - miniblob;
    %                     'Non-bridged blob removed!'
                    end
                end
                % Finalize miniblobs connected by bridge
                slice_connectedbybridge = bwlabel(logical(slice_connectedbybridge),4);
                n_slice_connectedbybridge = max(slice_connectedbybridge,[],'all');
                % Get information about each of these miniblobs
                miniblob_Zx_means = zeros(n_slice_connectedbybridge,1);
                miniblob_Zy_means = zeros(n_slice_connectedbybridge,1);
                miniblob_m_means  = zeros(n_slice_connectedbybridge,1);
                miniblob_mean_m   = zeros(n_slice_connectedbybridge,1);
                miniblob_directions = zeros(n_slice_connectedbybridge,1);
                miniblob_substantial= zeros(n_slice_connectedbybridge,1);
                miniblob_slanted  = zeros(n_slice_connectedbybridge,1);
                miniblob_rounded  = zeros(n_slice_connectedbybridge,1);
                for jj = 1:n_slice_connectedbybridge    % For each connected miniblob
                    miniblob = (slice_connectedbybridge==jj);
                    where_miniblob = find(miniblob);
                    miniblob_Zx_means(jj) = nanmean(slice_Zx(where_miniblob),'all');
                    miniblob_Zy_means(jj) = nanmean(slice_Zy(where_miniblob),'all');
                    miniblob_mean_m(jj)  = nanmean(sqrt(slice_Zx(where_miniblob).^2 + slice_Zy(where_miniblob).^2));     % Means of magnitude of slope. (e.g. Sphere treated as "slanted" with direction implied to be reliable.)
                    miniblob_m_means(jj)  = sqrt(nanmean(slice_Zx(where_miniblob),'all').^2 + nanmean(slice_Zy(where_miniblob),'all').^2);  % Magnitude of means of slope. (e.g. Sphere treated as "flat", so direction is unreliable.)
                    miniblob_directions(jj) = meanangle(-angle(slice_Zy(where_miniblob) - 1i*slice_Zx(where_miniblob)) * 180/pi,'all'); % Mean of angles of blob, CCW from North.
%                     miniblob_directions(jj) = -angle(miniblob_Zy_means(jj) - 1i*miniblob_Zx_means(jj)) * 180/pi;            % Angle from mean slopes of blob, CCW from North.
                    miniblob_substantial(jj)= sum(miniblob,'all') > (nansum(logical(slice_connectedbybridge),'all')*0.05); % If miniblob is >5% of total bridged region
                    
                    if miniblob_substantial(jj)
                        miniblob_slanted(jj)  = miniblob_m_means(jj) > slope_minslant;
                        miniblob_rounded(jj)  = miniblob_mean_m(jj) > 2.5*miniblob_m_means(jj);
                    else
                        miniblob_slanted(jj)  = miniblob_m_means(jj) > 5*slope_minslant;
                        miniblob_rounded(jj)  = miniblob_mean_m(jj) > 3*miniblob_m_means(jj);
                        % Small chunks have a less reliable slope, so they
                        % need to be steeper to be distinguished from the 
                        % rest of the slice.
                    end
                end
                anglediffs = @(angles) abs(angles - angles') - 360*(abs(angles - angles')>180);
                % FINAL CRITERIA:
                if (max(abs(anglediffs(miniblob_directions)),[],'all') > 70 ...
                 && ~ismember(false,miniblob_slanted)) ...  % If the slopes are distinctly rotated, and not because one is super flat...
                 || (ismember(true,miniblob_rounded)) ... % If one of the miniblobs is "rounded" (and therefore needs to be isolated)
                 || ((3*min(miniblob_m_means) < max(miniblob_m_means)) && ~ismember(false,miniblob_substantial))
                 % ^ ... or if one of them is flatter than the other (but not necessarily "~miniblob_slanted"), and both are of substantial size...
                    slice_bridge = slice_bridge + potential_bridge;  % ... then it's a "bridge"...
                    mask_split = mask_split | potential_bridge;
                end
            end
        end
        slice_split = slice .* ~slice_bridge;                   % ... and we'll cut any bridges out of `slice`!
    else
        slice_split = slice;
    end
end