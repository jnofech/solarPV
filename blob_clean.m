function [blobs] = blob_clean(binary,roof_mindiameter)
%BLOBS_CLEAN Given a binary matrix (1 indicating surface) and the smallest
%   possible roof size, returns labeled blobs with small structures
%   removed.
    binary = bwareaopen(logical(binary),ceil(roof_mindiameter^2),4);
    r = ceil(roof_mindiameter);         % Pixel "radius" of largest circle within smallest possible rooftop (More accurate, but restrictive)
%     r = ceil(sqrt(roof_mindiameter^2/pi));% Pixel "radius" of largest circle within smallest possible rooftop
    map_blobs = bwlabel(binary,4);      % Labeled blobs, including tiny structures.
    se = strel('disk',r);
    
    blobs       = zeros([size(binary)]);    % Final blob map (to be filled).
    for i = 1:max(map_blobs,[],'all')
        slice_full = (map_blobs==i);        % Binary array of each blob
        [rows, cols] = find(slice_full);
        row1 = min(rows); row2 = max(rows);
        col1 = min(cols); col2 = max(cols);

        slice = slice_full(row1:row2,col1:col2);        % Cropped array containing each blob
    
        % Delete small blobs ('Open' with circular SE fitting within smallest legally-allowed rooftop)
        slice_opened = imopen(slice,se);   % Any surface that exceeds bare minimum roof area.
%         if i==max(map_blobs,[],'all')
%             warning("SE radius is that of a circle considering panel length AND width. Maybe just the smaller one?");
%         end
        if sum(slice_opened,'all')==0
            slice = zeros(size(slice));
        end

        % If structure remains, add it to "blobs"
        % (Note: Adding instead of replacing allows blobs to be added
        % without zeros interfering with nearby blobs.)
        blobs(row1:row2,col1:col2) = blobs(row1:row2,col1:col2) + slice*i;
    end
    % Flatten and re-label blob map
    blobs = logical(blobs);
    blobs = bwareaopen(blobs,15);   % Cleans buggy "specks" that are left over
    blobs = bwlabel(blobs,4);
end

