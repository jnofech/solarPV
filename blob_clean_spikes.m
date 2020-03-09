function [blobs_bad] = blob_clean_spikes(isroof, Z)
%BLOB_CLEAN_SPIKES() identifies regions on rooftops where the height
%   may be changing very suddenly. These regions are not removed from
%   the blob, but will be ignored when calculating heights / slopes.
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
%   blobs_bad : double (2D array)
%       Labeled to each individual blob, indicating where heights suddenly
%       increase/decrease. These regions should ideally not be considered
%       when calculating roof heights/slopes.
    
    % ~~~ Initial Settings ~~~
    % Step size for "percentile vector" of heights. (Must be factor of 100.)
    percent_step = 0.5;        % ~0.5% is nice for detecting slopes
    binary = isroof;
    binary(isnan(binary)) = 0;

    map_blobs = bwlabel(isroof,4);          % Be aware that some "blobs" have literally 0 pixels.
    blobs_bad   = zeros(size(binary));    % 2D binary indicating bad regions of blobs (numbered).
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
        v_std = std(slice_Z,[],'all','omitnan');
        v_mean = nanmean(slice_Z,'all');
        vx = gradient(v,x)';

    %   % Conditions
        m_steep      = 10*nanmean(vx);  % Maximum slope of "percentile" vector (in m/%) at which surfaces are considered "steep".
        v_bar = 3*v_std;         % Max difference from mean height, at which surface is deemed "too high" (or too low)

        % Find percentiles (and heights) at which condition is met.
        % A percentile is selected if an entry on either side (NOT both sides)
        % meets condition.
        v_left  = circshift(v,-1); v_left(end) = 2*v_left(end-1) - v_left(end-2);
        v_right = circshift(v,1); v_right(1) = 2*v_right(2) - v_right(3);
        i_steep = (((v_left - v)/percent_step > m_steep) | ((v - v_right)/percent_step > m_steep))';
        i_toohigh = ((v > v_mean+v_bar) | (v < v_mean-v_bar))';

        % Finding specific height ranges where condition is met
        i_cond = i_toohigh | i_steep;            % 1 where condition on `v` is met.
        [1,diff(i_cond),1];                      % Nonzero where condition on `v` _changes_.
        idx = find([1,diff(i_cond),1]);          % Indices of the above. Doesn't indicate a condition ending until one space after.
        seq_lengths = idx(2:end)-idx(1:end-1);   % Step sizes over which condition is/isn't met. Will be looped over.
        j = 1;
        slice_highlight = zeros(size(slice));    % Binary array where condition is met.
        for step = seq_lengths
            if i_cond(j)==1 % & step
                heights = v( (j) : (j+step-1) );  % Examine the heights in this part of the sequence
                % Add to a slice-sized binary array wherever slice is within
                % `heights`.
                slice_highlight = slice_highlight + (slice_Z > min(heights)) & (slice_Z <= max(heights));
            end
            j = j+step;
        end
        % Use completed `slice_highlight` to add to `blobs_bad`.
%         blobs_bad(row1:row2,col1:col2) = blobs_bad(row1:row2,col1:col2) + slice_highlight*i;    % Old version: label matches blob
        blobs_bad(row1:row2,col1:col2) = blobs_bad(row1:row2,col1:col2) | slice_highlight;      % New version: binary
    end
end

