function roofs = blob_heights_split(isroof, Z)
%blob_heights_split() splits rooftops into multiple, smaller blobs if
%   multiple flat regions are detected.
%
%   Parameters:
%   -----------
%   isroof : logical (2D array)
%       Indicates any surface large enough to hold a tiny solar panel.
%   Z : double (2D array)
%       Height map of the building.
%
%   Returns:
%   --------
%   roofs : double (2D array)
%       `isroof`, but counting through the roofs.
    
    % ~~~ Initial Settings ~~~
    % Step size for "percentile vector" of heights. (Must be factor of 100.)
    percent_step = 2.5;        % 2.5% good for horizontal faces
    percent_step_split = 0.5;  % ~0.5% is nice for detecting slopes

    map_blobs = bwlabel(isroof,4);          % Be aware that some "blobs" have literally 0 pixels.
    blobs     = zeros([size(isroof)]);      % To be filled.
    n_blobs   = max(map_blobs,[],'all');
    
    for i = 1:n_blobs
        slice_full = (map_blobs==i);        % Binary array of each blob
        [rows, cols] = find(slice_full);
        row1 = min(rows); row2 = max(rows);
        col1 = min(cols); col2 = max(cols);

        slice = slice_full(row1:row2,col1:col2);        % Cropped binary array containing each blob
        slice_mask = double(slice);
        slice_mask(slice==0) = nan;                     % nan where roof isn't
        slice_Z = Z(row1:row2,col1:col2).*slice_mask;   % Heights of blob
        
        % Get "percentile vector" of heights
        x = percent_step:percent_step:100;
        v = prctile(slice_Z,x,'all');

        m_flat      = 0.17/25;   % Maximum slope of "percentile" vector (in m/%) at which surfaces are considered "flat".
        m_ultraflat = 0.05/25;   % Maximum slope of "percentile" vector (in m/%) at which surfaces are considered "steep".

        % Find percentiles (and heights) at which condition is met.
        % A percentile is selected if an entry on either side (NOT both sides)
        % meets condition.
        v_left  = circshift(v,-1); v_left(end) = 2*v_left(end-1) - v_left(end-2);
        v_right = circshift(v,1); v_right(1) = 2*v_right(2) - v_right(3);
        i_flat      = (((v_left - v)/percent_step < m_flat) | ((v - v_right)/percent_step < m_flat))';
        i_ultraflat = (((v_left - v)/percent_step < m_ultraflat) | ((v - v_right)/percent_step < m_ultraflat))';
%         i_steep = (((v_left - v)/percent_step > m_steep) | ((v - v_right)/percent_step > m_steep))';

        % Finding specific height ranges where condition is met
        i_cond = i_flat;                         % 1 where condition on `v` is met.
        [1,diff(i_cond),1];                      % Nonzero where condition on `v` _changes_.
        idx = find([1,diff(i_cond),1]);          % Indices of the above. Doesn't indicate a condition ending until one space after.
        seq_lengths = idx(2:end)-idx(1:end-1);   % Step sizes over which condition is/isn't met. Will be looped over.
        j = 1;
        flat_count = 0;
        for step = seq_lengths
            if i_cond(j)==1 % & step
                heights = v( (j) : (j+step-1) );  % Examine the heights in this part of the sequence
                % If this part of the sequence spans ___% of the percentile
                % array, then it's deemed a "flat surface". 
                if (step-1)*percent_step >= 5.00
                    flat_count = flat_count + 1;
                    if flat_count==1
                        x_inbetween_min = x(j+step-1);
                        v_inbetween_min = v(j+step-1);
                    end
                    x_inbetween_max = x(j);
                    v_inbetween_max = v(j);
                end
            end
            j = j+step;
        end
        % If following conditions are met, split the slice!
        if flat_count >= 2 ...      % Multiple flat regions...
        & (v_inbetween_max-v_inbetween_min > 1.0) ...  % ... with a minimum height difference...
        & ((v_inbetween_max-v_inbetween_min)/(x_inbetween_max-x_inbetween_min) > 2*m_flat)  % ... and the heights in between are not too flat... 
            disp("Slice "+i+": Multiple flat surfaces ("+flat_count+") detected.");
            % Generate new "percentile vector", with much smaller steps this time.
            x_split = percent_step_split:percent_step_split:100;
            v_split = prctile(slice_Z,x_split,'all');
            vx_split = gradient(v_split,x_split)';
%             [vx_pks,x_pks] = findpeaks(vx_split,x_split);   % Requires Signal Processing Toolkit?
            is_pks = islocalmax(vx_split);
            vx_pks = vx_split(is_pks); x_pks = x_split(is_pks);
            is_inbetween = (x_pks > x_inbetween_min) .* (x_pks < x_inbetween_max);
            vx_peak = max(vx_pks.*is_inbetween);  % Peak "steepness" of percentile vector, between flat surfaces.
            where_peak = find(vx_pks == vx_peak);
            % Mark on "v" where slope is within some value of the peak slope
            % (i.e. find the range of heights to be treated as "bad").
            % Delete all heights in this range from slice.
            n_split_miniblobs = 1;      % Cut parts from `slice` until >=2 miniblobs appear.
            vx_width      = 0.050;
            vx_width_step = 0.025;
            while (n_split_miniblobs == 1) && (vx_width < vx_peak*2)
                where_nearpeak = find(  vx_split > (vx_peak-vx_width) ...
                                      & vx_split < (vx_peak+vx_width) ... % <- redundant; slope can't get higher than peak...
                                      & x_split  > x_inbetween_min ...
                                      & x_split  < x_inbetween_max  );
                v_nearpeak = v_split(where_nearpeak);
                slice_remove = (slice_Z > min(v_nearpeak) & slice_Z < max(v_nearpeak)); % Chunks of slice to be removed

                slice_split = slice .* ~slice_remove;    % Binary of slice, with certain heights removed
                slice_split = bwareaopen(slice_split, round(nansum(slice,'all')*0.025), 4); % Remove small chunks.
                slice_split_miniblobs = bwlabel(slice_split,4);
                n_split_miniblobs = max(slice_split_miniblobs,[],'all');
                vx_width = vx_width+vx_width_step;
            end
            
            % If a split has occurred, then make sure only the "connecting"
            % area is removed!
            if n_split_miniblobs >= 2
                % Identify which parts of `slice_remove` act as "bridges" between miniblobs, and only consider those
                slice_bridge = zeros(size(slice_split)); % Miniblobs that join pieces of `slice_split` together
                slice_remove_miniblobs = bwlabel(slice_remove,4);
                n_remove_miniblobs = max(slice_remove_miniblobs,[],'all');
                for ii = 1:n_remove_miniblobs   % Check each `slice_remove` blob to see whether it's a "bridge"
                    remove_miniblob = (slice_remove_miniblobs == ii);
                    n_miniblobs_bridged = max(bwlabel((slice_split + remove_miniblob),4),[],'all');
                    if n_miniblobs_bridged < n_split_miniblobs    % If this miniblob "connects" two parts of `slice_split`...
                        slice_bridge = slice_bridge + remove_miniblob;  % ... then it's a "bridge"...
                    end
                end
                slice_split = slice .* ~slice_bridge;                   % ... and we'll cut any bridges out of `slice`!
                % Save changes into slice!
                slice = slice_split;
                slice_mask = double(slice);
                slice_mask(slice==0) = nan;         % nan where roof isn't
            end
        end
        
%         % Plot specific blob
%         if i==6
%             fontsize = 14;
%             figure;
%             width=1600;height=900;
%             set(gcf,'position',[160,60,width,height]);
%             subplot(1,2,1); imagesc(slice_Z); axis image; colorbar; title("Heights of slice (m)");
%             ax = gca; ax.FontSize = fontsize;
%             subplot(1,2,2); plot(x,v,'blue'); hold on; ... % plot(x,v_left,'red'); plot(x,v_right,'green'); ...
%                 x_cond = x(find(i_cond)); y_range = ylim(); x_nearpeak = x_split(where_nearpeak);  ...
%                 plot([x_cond; x_cond], repmat(ylim',1,size(x_cond,2)), '-k'); ...
%                 plot([x_nearpeak; x_nearpeak], repmat(ylim',1,size(x_nearpeak,2)), '-r'); ...
%                 ylim(y_range); hold off;
%                 title("Percentile analysis"); ylabel("Height percentile (m)"); xlabel("Percent rank");
%                 ax = gca; ax.FontSize = fontsize;
% %             subplot(2,2,4); plot(x,i_cond); xlabel("Percentile"); title("Is flat?");
% %             ax = gca; ax.FontSize = fontsize;
%         end
        
        % Add to map of blobs
        blobs(row1:row2,col1:col2) = blobs(row1:row2,col1:col2) + slice*1;
    end
    blobs = bwareaopen(blobs,2,4);
    blobs = bwlabel(logical(blobs),4);
    roofs = blobs;
end

