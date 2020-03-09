function [roofs_modified,roof_isbad_modified,blob_info] = blob_classify_noslopesplit(roofs,roof_isbad,Z,dZdX,dZdY,pix2m,roof_mindiameter,loop)
%BLOB_CLASSIFY Classifies each blob as "flat" (1), "slanted" (2), or
%   "irregular" (3); and returns height, parapet, and/or slope information
%   accordingly.
% 
%   Parameters:
%   -----------
%   roofs : logical (2D array)
%       Labeled blobs.
%   roofs_bad : double (2D array)
%       "Unreliable" regions in each labeled blob.
%   Z : double (2D array)
%       Height map of the building.
%   dZdX : double (2D array)
%       Slopes (East) of the building.
%   dZdY : double (2D array)
%       Slopes (South) of the building.
%   pix2m : double
%       (TEMP!!) Metres per pixel
%   roof_minsize_pix : double
%       Minimum area of rooftop, in pixels.
%   loop : double (integer)
%       Loop counter. Will stop classifying rooftops after the 5th loop.
%
%   Returns:
%   --------
%   roofs_modified : logical (2D array)
%       Labeled blobs, after modifications.
%   roof_isbad_modified : logical (2D array)
%       Binary indicating unreliable "spikes" on "roofs_modified".
%   blob_info : table
%       `roof`          : Blob number of rooftop
%       `type`          : Blob classification
%       `height_roof`   : Height at highest point of roof
%       `height_parapet`: Approximate parapet height (if any)
%       `slope_x`       : Slope in X-direction (east?) (if any)
%       `slope_y`       : Slope in Y-direction (south?) (if any)
    loop_max = 4;       % Any iteration reaching this point will be completely "irregular".
    n_blobs = max(roofs,[],'all');
    
    % Set up "irregular" blobs
    roofs_irr = roofs;          % (will be) `roofs` for irregular blobs
    roofs_copy = roofs;         % (will be) `roofs` with irregulars deleted
    
    % Set up table
    roof = (1:n_blobs)';
    type = strings(n_blobs,1);
    height_roof = zeros(n_blobs,1);
    height_parapet = zeros(n_blobs,1);
    slope_x = zeros(n_blobs,1);
    slope_y = zeros(n_blobs,1);
    slope   = zeros(n_blobs,1);
    direction = nan(n_blobs,1);
    area = zeros(n_blobs,1);
    binary = cell(n_blobs,1);
    
    for i = 1:n_blobs
        slice_full = (roofs==i).*~(roof_isbad);         % Binary array of each blob, EXCLUDING "bad" regions
        slice_full_full = (roofs==i);                   % Binary array of each blob, INCLUDING "bad" regions (for dimensions/area considerations)
        [rows, cols] = find(slice_full);
        row1 = min(rows); row2 = max(rows);
        col1 = min(cols); col2 = max(cols);

        slice   = slice_full(row1:row2,col1:col2);      % Cropped array containing each blob, EXCLUDING "bad" regions
        slice_mask = double(slice);
        slice_mask(slice==0) = nan;                     % nan where roof isn't
        slice_Z = Z(row1:row2,col1:col2).*slice_mask;   % Heights of blob

        [isflat,height_roof(i),height_parapet(i)] = slice_classify_isflat(slice_Z);

        if isflat && (loop~=loop_max)
            % It's flat and horizontal.
            type(i) = "flat";
            slope_x(i) = 0;
            slope_y(i) = 0;
            if height_parapet(i) < 0.20
                height_parapet(i) = 0;
            end
            roofs_irr      = roofs_irr - i*(roofs==i);
        else
            height_parapet(i) = 0;
            [isslanted,height_roof(i),slope_x(i),slope_y(i)] = slice_classify_isslanted(slice_Z);
            if isslanted && (loop~=loop_max)
                % It's flat and slanted.
                type(i) = "slanted";
                slope_m  = sqrt(slope_x(i)^2 + slope_y(i)^2);
                slope(i) = atand(slope_m);
                direction(i) = -angle(slope_y(i) + 1i*slope_x(i)) * 180/pi;  % Angle of building, CCW from North.
                roofs_irr      = roofs_irr - i*(roofs==i);
            else
                % It's irregular.
                type(i) = "delete";
                height_parapet(i) = nan;
                slope(i) = nan;
                direction(i) = nan;
                roofs_copy = roofs_copy - i*(roofs==i);   % Irregular roofs will be re-added later.
                                                          % (!!!) roof_isbad will be
                                                          % masked accordingly based
                                                          % on final result.
            end
        end
            
        % TEMPORARY CHANGES FOR NIMA
        area(i) = sum(slice_full_full,'all') * pix2m^2;
        if isfinite(slope(i))
            area(i) = area(i) / cosd(slope(i));
        end
        binary{i} = slice_full_full(row1:row2,col1:col2);
    end
    
    % Remove "deleted" roofs from table.
    disp("Number of `irregular` roofs: "+string(sum(type=="delete"))+" out of "+size(type,1));
%     if sum(type=="delete")==size(type,1)
    if true
        % If they're ALL `deleted`, or if it's looped too many times:
        type(type=="delete") = "irregular";
        roofs_modified = roofs;
        roof_isbad_modified = roof_isbad;   % No modifications.
        blob_info = table(roof,type,height_roof,height_parapet,slope,direction,binary,area);
        return;
%     else
%         % If there are still not-`irregular` roofs:
%         roofs = roofs_copy;
%         roof(type=="delete") = [];
%         height_roof(type=="delete") = [];
%         height_parapet(type=="delete") = [];
%         slope(type=="delete") = [];
%         direction(type=="delete") = [];
%         binary(type=="delete") = [];
%         area(type=="delete") = [];
%         type(type=="delete") = [];
%         % Re-order the current "good" roofs
%         roofs = bwlabel(logical(roofs),4);
%         n_blobs = max(roofs,[],'all');
%         roof = (1:n_blobs)';
%         if size(roof,1)~=n_blobs
%             error('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA');
%         end
    end
    
    % Fix irregular roofs
    roofs_irr = bwlabel(logical(roofs_irr),4);
    isroof_irr = regional_minmax_split(Z,dZdX,dZdY,roofs_irr,pix2m);
    isroof_irr = blob_clean(isroof_irr,roof_mindiameter);   % Cleaning
    roof_isbad_irr = roof_isbad .* isroof_irr;
    roofs_irr = bwlabel(isroof_irr,4);
    
    % Classify the improved irregular roofs!
    % (When the `blob_classify` loop is broken, ALL of these will be
    % type `delete`. Final result should be mixture of good AND bad roofs.)
    [roofs_irr,~,info_irr] = blob_classify(roofs_irr,roof_isbad_irr,Z,dZdX,dZdY,pix2m,roof_mindiameter,loop+1);
    
    % Add the results on top of previous table!
    roof            = cat(1,roof,info_irr.roof+n_blobs);
    height_roof     = cat(1,height_roof,info_irr.height_roof);
    height_parapet  = cat(1,height_parapet,info_irr.height_parapet);
    slope           = cat(1,slope,info_irr.slope);
    direction       = cat(1,direction,info_irr.direction);
    binary          = cat(1,binary,info_irr.binary);
    area            = cat(1,area,info_irr.area);
    type            = cat(1,type,info_irr.type);
    
    % Finally, add the improved roofs_irr and roofs_isbad_irr to the rest,
    % and set up the table!
    roofs_modified = roofs + (roofs_irr+logical(roofs_irr)*n_blobs);  % Starts counting irregulars at end
    roof_isbad_modified = roof_isbad .* logical(roofs_modified);
    blob_info = table(roof,type,height_roof,height_parapet,slope,direction,binary,area);
end

function [isflat,height_roof,dzdx,dzdy] = slice_classify_isslanted(slice_Z)
% Checks if slice is flat, and returns plane-fitted height and slopes.
% Note that this includes horizontal surfaces as well, although it is
% currently better to just use this for slanted surfaces.
    u = symunit;

    % Create X, Y arrays
    indices = find(slice_Z==nanmax(slice_Z,[],'all'));
    [ymax,xmax] = size(slice_Z);
    [rows,cols] = ind2sub(size(slice_Z),indices);
    row_cen = round(mean(rows));    col_cen = round(mean(cols));
    x = 1:xmax; y = 1:ymax; [X,Y] = meshgrid((x-col_cen)*val(rewrite(u.pix,u.m)),(y-row_cen)*val(rewrite(u.pix,u.m)));
    
    % Convert X, Y, Z to 1D and remove NaNs
    x = X(:);   y = Y(:);   z = slice_Z(:);
    x = x(~isnan(z)); y = y(~isnan(z)); z = z(~isnan(z));
    I = ones(size(x));
    A = [I x y];
    B = z;
    coeff = A\B;                % Least-squares fit of a plane to Z.
    height_roof = coeff(1);     % Peak of rooftop
    dzdx        = coeff(2); 
    dzdy        = coeff(3);
    Zfit = coeff(1) + coeff(2)*X + coeff(3)*Y;
    Zdiff = slice_Z - Zfit;
    stdev = std(Zdiff,[],'all','omitnan');
    
%     warning("stdev threshold for plane fitting arbitrarily set to 0.09. Need to test on sloped surfaces and irregular surfaces!");
%     if stdev < 0.20
    if stdev < 0.09
        isflat = true;
    else
        % Plane fit didn't go well
        isflat = false;
        dzdx = nan;
        dzdy = nan;
    end
end

function [isflat,height_roof,height_parapet] = slice_classify_isflat(slice_Z)
% Checks if slice is flat, and returns information accordingly.
    % ~~~ Initial Settings ~~~
    % Step size for "percentile vector" of heights. (Must be factor of 100.)
    percent_step = 2.5;        % 2.5% good for horizontal faces
    ultraflat_weight = 3;      % Weight multiplier for `ultraflat` regions, when averaging height
    
    % Get "percentile vector" of heights
    x = percent_step:percent_step:100;
    v = prctile(slice_Z,x,'all');  % v(end) = nanmax(slice_Z,[],'all');
    m_flat      = 0.17/25;   % Maximum slope of "percentile" vector (in m/%) at which surfaces are considered "flat".
    m_ultraflat = 0.05/25;   % Maximum slope of "percentile" vector (in m/%) at which surfaces are considered "steep".
    
    % Find percentiles (and heights) at which condition is met.
    % A percentile is selected if an entry on either side (NOT both sides)
    % meets condition.
    v_left  = circshift(v,-1); v_left(end) = 2*v_left(end-1) - v_left(end-2);
    v_right = circshift(v,1); v_right(1) = 2*v_right(2) - v_right(3);
    i_flat      = (((v_left - v)/percent_step < m_flat) | ((v - v_right)/percent_step < m_flat))';
    i_ultraflat = (((v_left - v)/percent_step < m_ultraflat) | ((v - v_right)/percent_step < m_ultraflat))';
    
    % Determine whether slice is flat
    i_cond = i_flat;                         % 1 where condition on `v` is met.
    isflat = false;
    height_roof = nan;
    height_parapet = nan;
%     if sum(i_cond)*percent_step > 50.00      % <-- determining factor for whether slice `isflat`
    if sum(i_cond)*percent_step > 50.00 ...
    && sum(i_ultraflat)*percent_step > 10.00 % <-- determining factor for whether slice `isflat`
        isflat = true;
        [1,diff(i_cond),1];                      % Nonzero where condition on `v` _changes_.
        idx = find([1,diff(i_cond),1]);          % Indices of the above. Doesn't indicate a condition ending until one space after.
        seq_lengths = idx(2:end)-idx(1:end-1);   % Step sizes over which condition is/isn't met. Will be looped over.
        j = 1;
        flat_count = 0;
        for step = seq_lengths
            if i_cond(j)==1 % If this sequence meets `i_cond`
                heights = v( (j) : (j+step-1) );  % Examine the heights in this part of the sequence
                j_ultraflat = i_ultraflat( (j) : (j+step-1) )';  % Indicates "ultraflat" regions in this sequence
                % If this part of the sequence spans ___% of the percentile
                % array, then it's deemed a "flat surface". 
                if (step-1)*percent_step >= 5.00
                    flat_count = flat_count + 1;  % Number of flat regions that have been found in this slice so far
                    if flat_count==1
                        height_lower = sum( (heights .* ~j_ultraflat) + ultraflat_weight*(heights .* j_ultraflat)) ...
                                       / sum( nnz(heights .* ~j_ultraflat) + ultraflat_weight*nnz(heights .* j_ultraflat));
                        % ^ Ultraflat region has greater weighting
                    end
                    height_upper = sum( (heights .* ~j_ultraflat) + ultraflat_weight*(heights .* j_ultraflat)) ...
                                   / sum( nnz(heights .* ~j_ultraflat) + ultraflat_weight*nnz(heights .* j_ultraflat));
                    % ^ Ultraflat region has greater weighting
                end
            end
            j = j+step;
        end
        
        % Parapet warning!
        if flat_count >= 2
%             height = height of the bottom flat region (OR the larger one, if the bottom one is too tiny)
%             parapet height = approximately (height of highest region - height of lower region)
            height_roof = height_lower;
            height_parapet = height_upper - height_roof;
        else
%             height = mean of all "flat" heights
%             parapet height = <-- Get rid of everything below the highest "flat" height. Take ~90th percentile of everything above it!
            height_roof = height_lower;
            slice_Z_parapet = slice_Z;
            slice_Z_parapet(slice_Z_parapet < max(heights)) = nan;  % Note: `heights` = range of heights for highest rooftop region
            height_parapet = prctile(slice_Z_parapet,90,'all') - height_roof;
        end
    end
end