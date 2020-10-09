function [roofs_modified,roofs_modified_p,roof_isbad_modified,blob_info] = blob_classify(roofs,roof_isbad,Z,dZdX,dZdY,pix2m,roof_mindiameter,loop)
%BLOB_CLASSIFY Classifies each blob as "flat" (1), "slanted" (2), or
%   "irregular" (3); and returns height, parapet, and/or slope information
%   accordingly.
% 
%   Parameters:
%   -----------
%   roofs : double (2D array)
%       Labeled blobs.
%   roofs_bad : logical (2D array)
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
%       `coords`        : Coordinates (row1, col1) of first point in binary

    loop_max = 6;       % Any iteration reaching this point will be completely "irregular".
%     loop_max = 1;       % Any iteration reaching this point will be completely "irregular".
    n_blobs = max(roofs,[],'all');
    
    % Set up roof-tracking blobs
    roofs_irreg = zeros(size(roofs));       % (will be) `roofs` for irregular blobs
    roofs_good   = zeros(size(roofs));      % (will be) `roofs` with irregulars deleted; parapets excluded
    roofs_good_p = zeros(size(roofs));      % (will be) `roofs` with irregulars deleted; parapets included
    isroof_irr_dontsplit = zeros(size(roofs));  % (will be) `logical(roofs_irreg)`, but for roofs that I don't want split
    isroof_toosmall = zeros(size(roofs));   % (will be) `logical(roofs)` for flat/slanted blobs 
                                            % that are too small to fit a panel
    
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
    toosmall = false(n_blobs,1);
    coords = cell(n_blobs,1);
    binary = cell(n_blobs,1);
    binary_parapet = cell(n_blobs,1);
    
    for i = 1:n_blobs
        slice_full     = (roofs==i).*~(roof_isbad); % Binary array of each blob, EXCLUDING "bad" regions
        slice_full_all = (roofs==i);                % Binary array of each blob, INCLUDING "bad" regions (for dimensions/area considerations)
        
        % Crop main maps!
        edge_buffer = 10;
        [slice_all,slice,slice_Z_unmasked] = crop(slice_full_all,slice_full,Z,edge_buffer);
            % Track map position (for "repopulating" `roofs_good` with modified roofs)
            [rows, cols] = find(slice_full_all);
            row1 = min(rows)-edge_buffer; row2 = max(rows)+edge_buffer;
            col1 = min(cols)-edge_buffer; col2 = max(cols)+edge_buffer;
            row1 = max(row1,1);     row2 = min(row2,size(roofs,1));
            col1 = max(col1,1);     col2 = min(col2,size(roofs,2));

        % Create masks for Z, so that background isn't caught in height analysis
        slice_mask = double(slice);
        slice_mask(slice==0)         = nan;     % nan where roof isn't
        slice_mask_all = double(slice_all);
        slice_mask_all(slice_all==0) = nan;     % nan where roof isn't

        slice_Z     = slice_Z_unmasked.*slice_mask;     % Heights, EXCLUDING "bad" regions
        slice_Z_all = slice_Z_unmasked.*slice_mask_all; % Heights, INCLUDING "bad" regions
        
        % Get classification!
        if loop~=loop_max
%             [type(i),height_roof(i),height_parapet(i),slope(i),direction(i)] = slice_classify_flatorslanted(slice,slice_all,slice_Z,slice_Z_all);
            [type(i),height_roof(i),height_parapet(i),slope(i),direction(i),slice_noparapet,slice_parapet] = slice_classify_flatorslanted(slice_all,slice_Z,slice_Z_all,pix2m);
        else
            % Out of loops; it's automatically irregular.
            type(i) = "delete";
            height_roof(i) = max(slice_Z,[],'all');
            height_parapet(i) = nan;
            slope(i) = nan;
            direction(i) = nan;
            
            % Define `slice_noparapet`, just in case:
            slice_noparapet = slice_all;
        end
        % Keep track of good/bad roofs.
        if type(i)=="flat" || type(i)=="slanted"
            % Check if it can fit a panel with a walkway around it.
            se = strel('disk',ceil(roof_mindiameter * cosd(slope(i))));
            slice_test = imerode(slice,se);
            % Check size
            roofs_good(row1:row2,col1:col2)   = roofs_good(row1:row2,col1:col2) + i*slice_noparapet;                    % No parapet
            roofs_good_p(row1:row2,col1:col2) = roofs_good_p(row1:row2,col1:col2) + i*(slice_noparapet+slice_parapet);  % Yes, parapet
            if sum(slice_test,'all')>0
                % It's good!
            else
                % It would be good, but is too small.
                toosmall(i) = true;
%                 isroof_toosmall = isroof_toosmall + roofs==i;
            end
            
%                             % Check size (old)
%                             if sum(slice_test,'all')>=0
%                                 % It's good!
%                 %                 roofs_good = roofs_good + i*slice_full_all; % Original roof
%                 %                 roofs_good(row1:row2,col1:col2) = roofs_good(row1:row2,col1:col2) + i*slice_all; % Original roof
%                                 roofs_good(row1:row2,col1:col2)   = roofs_good(row1:row2,col1:col2) + i*slice_noparapet;                    % No parapet
%                                 roofs_good_p(row1:row2,col1:col2) = roofs_good_p(row1:row2,col1:col2) + i*(slice_noparapet+slice_parapet);  % Yes, parapet
%                             else
%                                 % IRREGULAR, due to being too small
%                 %                 disp("REMOVED! Rooftop is too small.");
%                                 type(i) = "delete";
%                                 isroof_toosmall = isroof_toosmall + roofs==i;
%                             end
            
            % Check if a flat roof has split into multiple pieces due to
            % parapet detection.
            n_flats = max(bwlabel(slice_noparapet,4),[],'all');
            if n_flats>1
                disp("FLAT ROOF #"+i+" HAS SPLIT INTO MULTIPLE PIECES DUE TO PARAPET DETECTION");
                type(i) = "delete";
                % Remove them from the "good rooftop" tracking map
                roofs_good(row1:row2,col1:col2)   = roofs_good(row1:row2,col1:col2) - i*slice_noparapet;                    % No parapet
                roofs_good_p(row1:row2,col1:col2) = roofs_good_p(row1:row2,col1:col2) - i*(slice_noparapet+slice_parapet);  % Yes, parapet
                isroof_irr_dontsplit(row1:row2,col1:col2) = isroof_irr_dontsplit(row1:row2,col1:col2) + slice_noparapet;
            end
        else
            % IRREGULAR, due to actual classification
            type(i) = "delete";
            roofs_irreg = roofs_irreg + i*slice_full_all;
        end
            
        % (!) ADDING STATS TO TABLE
        area(i) = sum(slice_noparapet,'all') * pix2m^2;
%         area(i) = sum(slice_full_all,'all') * pix2m^2;
%         if isfinite(slope(i))
%             area(i) = area(i) / cosd(slope(i));
%         end
        % Track map position (Coords include parapets, but parapets are ignored in polygon approximation!)
            [rows, cols] = find((roofs_good_p==i) + (roofs_irreg==i));
            row1 = min(rows); row2 = max(rows);
            col1 = min(cols); col2 = max(cols);
        coords{i} = [row1,col1];
        binary{i}         = (roofs_good(row1:row2,col1:col2)==i) + (roofs_irreg(row1:row2,col1:col2)==i);    % a.k.a. `binary_noparapet`
        binary_parapet{i} = (roofs_good_p(row1:row2,col1:col2)==i) + (roofs_irreg(row1:row2,col1:col2)==i);    % a.k.a. `binary_yesparapet`
    end
    
                    %     figure(1); imagesc(logical(roofs_good) + 5*logical(roofs_irreg)); axis image; colorbar
                    %     figure(2); imagesc(roofs_good + roofs_irreg); axis image; colorbar
    
    % Remove "deleted" roofs from table, if splitting is to be attempted.
    disp("Number of `irregular` roofs: "+string(sum(type=="delete"))+" out of "+size(type,1));
    if sum(type=="delete")==size(type,1)
        % If they're ALL `deleted`, or if it's looped too many times:
        type(type=="delete") = "irregular";
        roofs_modified   = roofs;
        roofs_modified_p = roofs;
        roof_isbad_modified = roof_isbad;   % No modifications.
%         imagesc(roofs); axis image; colorbar; title('BAD roofs');
        blob_info = table(roof,type,height_roof,height_parapet,slope,direction,coords,binary,binary_parapet,area,toosmall);
        return;
    else
        % If there are still not-`irregular` roofs:
        roofs   = roofs_good;
        roofs_p = roofs_good_p;
        roof(type=="delete") = [];
        height_roof(type=="delete") = [];
        height_parapet(type=="delete") = [];
        slope(type=="delete") = [];
        direction(type=="delete") = [];
        coords(type=="delete") = [];
        binary(type=="delete") = [];
        binary_parapet(type=="delete") = [];
        area(type=="delete") = [];
        toosmall(type=="delete") = [];
        type(type=="delete") = [];
        % Re-order the current "good" roofs
        roofs_reordered = bwlabel(logical(roofs),4);
        roofs_p_reordered = zeros(size(roofs)); % (Will be) roofs with parapets, reordered to correspond with `roofs_reordered`.
        n_roofs_reordered = max(roofs_reordered,[],'all');
        for i = 1:n_roofs_reordered
            loc = find(roofs_reordered==i);
            loc = loc(1);   % Location of some point inside roof `i`.
            old_index = roofs(loc);         % Old index of roof `i`.
            roofs_p_reordered = roofs_p_reordered + i*(roofs_p==old_index);
        end
        roofs   = roofs_reordered;
        roofs_p = roofs_p_reordered;
%         roofs = bwlabel(logical(roofs),4);
        n_blobs = max(roofs,[],'all');
        roof = (1:n_blobs)';
        if size(roof,1)~=n_blobs
            error('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA');
        end
    end
    
    % Fix irregular roofs
    roofs_irreg = bwlabel(logical(roofs_irreg),4);
    isroof_irr = regional_minmax_split(Z,dZdX,dZdY,roofs_irreg,roof_isbad,pix2m); % Slope-based splitting
%     isroof_irr = isroof_irr + isroof_toosmall;  % Adding "too-small" roofs back in
    isroof_irr = isroof_irr + isroof_irr_dontsplit;  % Adding more irregular roofs, but without attempting splitting on them
    isroof_irr = blob_clean(isroof_irr,roof_mindiameter/2);   % Cleaning
    roof_isbad_irr = roof_isbad .* isroof_irr;
    roofs_irreg = bwlabel(isroof_irr,4);
    disp("  After splitting roofs, there are now "+max(roofs_irreg,[],'all')+" irregular rooftops.");
    
    % Classify the improved irregular roofs!
    % (When the `blob_classify` loop is broken, ALL of these will be
    % type `delete`. Final result should be mixture of good AND bad roofs.)
    [roofs_irreg,~,~,info_irr] = blob_classify(roofs_irreg,roof_isbad_irr,Z,dZdX,dZdY,pix2m,roof_mindiameter,loop+1);
    
    % Add the results on top of previous table!
    roof            = cat(1,roof,info_irr.roof+n_blobs);
    height_roof     = cat(1,height_roof,info_irr.height_roof);
    height_parapet  = cat(1,height_parapet,info_irr.height_parapet);
    slope           = cat(1,slope,info_irr.slope);
    direction       = cat(1,direction,info_irr.direction);
    coords          = cat(1,coords,info_irr.coords);
    binary          = cat(1,binary,info_irr.binary);
    binary_parapet  = cat(1,binary_parapet,info_irr.binary_parapet);
    area            = cat(1,area,info_irr.area);
    toosmall        = cat(1,toosmall,info_irr.toosmall);
    type            = cat(1,type,info_irr.type);
    
    % Finally, add the improved roofs_irr and roofs_isbad_irr to the rest,
    % and set up the table!
    roofs_modified   = roofs   + (roofs_irreg+logical(roofs_irreg)*n_blobs);  % Starts counting irregulars at end
    roofs_modified_p = roofs_p + (roofs_irreg+logical(roofs_irreg)*n_blobs);  % Starts counting irregulars at end
    roof_isbad_modified = roof_isbad .* logical(roofs_modified);
    blob_info = table(roof,type,height_roof,height_parapet,slope,direction,coords,binary,binary_parapet,area,toosmall);
end

function [classification,height_roof,height_parapet,slope,direction,slice_all_noparapet,slice_parapet] = slice_classify_flatorslanted(slice_all,slice_Z,slice_Z_all,pix2m)
% Check if slice is flat.
% If slice has significant "flat region":
%   Classified as 'flat'. Slope=0; direction=NaN.
% else
%   Check if slice is slanted.
%   if slice is slanted:
%       if slope is <0.02 (i.e. minimum slope for 
%         Classified as 'flat'. Slope=0; direction=NaN.
%       else
%         Classified as 'slanted'. 
%         Calculate slope+direction. Parapet height = NaN.
%   else
%       Classified as 'irregular'.
%       Slope, direction, parapet height = NaN.
    [isflat,height_roof] = slice_classify_isflat(slice_Z);
    if isflat
        classification = "flat";
        slope = 0;
        direction = NaN;
    else
        [isslanted,height_roof,slope_x,slope_y] = slice_classify_isslanted(slice_Z);
        if isslanted
            slope_m  = sqrt(slope_x^2 + slope_y^2);
            if slope_m < 0.02
                classification = "flat";
                slope = 0;
                direction = NaN;
            else
                classification = "slanted";
                slope = atand(slope_m);
                direction = -angle(slope_y + 1i*slope_x) * 180/pi;  % Angle of building, CCW from North. -180 to +180.
            end
        else
            classification = "irregular";
            slope = nan;
            direction = nan;
        end
    end
    
    % Parapet stuff!
    if isflat
        [slice_all_noparapet,slice_parapet,height_parapet] = slice_findparapet(slice_all,height_roof,slice_Z_all,pix2m);
    else
        slice_all_noparapet = slice_all;
        slice_parapet = zeros(size(slice_all));
        height_parapet = NaN;
    end 
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

function [isflat,height_roof] = slice_classify_isflat(slice_Z)
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
%     disp("Percent flat: "+sum(i_cond)*percent_step)
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
        
        % (OLD) Parapet warning!
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
    
%     % Plot for a specific rooftop
%     if height_roof >= 35.4475 && height_roof <= 35.4477 % Top roof of CCIS
%     if height_roof >= 13.283 && height_roof <= 13.285 % Pembina?
%         figure;
%         width=1600;height=900; fontsize = 14;
%         set(gcf,'position',[160,60,width,height]);
%         subplot(1,2,1); imagesc(slice_Z); axis image; colorbar;
%         ax = gca; ax.FontSize = fontsize;
%             if isflat==true
%                 title("FLAT SURFACE;"+newline ...
%                       +"roof height  = "+height_roof+", "+newline ...
%                       +"para. height = "+height_parapet);%+" ("+(height_parapet+height_roof)+")");
%             else
%                 title(newline ...
%                       +"roof height  = "+height_roof+", "+newline ...
%                       +"para. height = "+height_parapet);%+" ("+(height_parapet+height_roof)+")");
%             end
%         subplot(1,2,2); plot(x,v,'blue'); hold on; ... %plot(x,v_left,'red'); plot(x,v_right,'green'); ...
%     %     i_cond = i_ultraflat;                         % 1 where condition on `v` is met.
%             x_cond = x(find(i_cond)); y_range = ylim(); ...
%             plot([x_cond; x_cond], repmat(ylim',1,size(x_cond,2)), '-k'); ...
%             ylim(y_range); hold off;
%             title("Percentile analysis"); ylabel("Height percentile (m)"); xlabel("Percent rank");
%             ax = gca; ax.FontSize = fontsize;
% %         subplot(2,2,4); plot(x,i_cond); xlabel("Percentile"); title("Meets condition?");
% %         ax = gca; ax.FontSize = fontsize;
%     end
end

function [slice_all_noparapet,slice_parapet,height_parapet] = slice_findparapet(slice_all,height_roof,slice_Z_all,pix2m)
    % Fill up holes and cracks of roof
    dilation_size = 0.8;    % in m
    dilation_size_pix = ceil(dilation_size / pix2m);
    se = strel('disk',dilation_size_pix);
    slice_fat = imdilate(slice_all,se);
    slice_fat = imfill(slice_fat,'holes');
    slice_all_filled = imerode(slice_fat,se) | slice_all;

    % Parapet assumptions (from
    % http://www.infrastructure.alberta.ca/Content/docType486/Production/redbook-14thedition.pdf):
    %   Width (thickness): at least 0.38m. I'll go with a max of ~1m.
    %   Height: minimum 30in, from Google (76.2cm). I'll use min of 60cm, in
    %   case Google models aren't "extreme" enough.
    para_maxwidth = 1.0;    % Maximum parapet width, in m.
    para_maxwidth = ceil(para_maxwidth / pix2m);
    para_minheight = 0.60;  % Minimum parapet height, in m.

    % Separate "outer roof" from "inner roof". Note that inner roof may
    % disappear entirely if it's <~`2*para_width` wide.
    slice_inner = imerode(slice_all_filled,strel('disk',para_maxwidth));
    slice_outer = slice_all_filled - slice_inner;

    % Find parapet! And then smooth it out, by dilating->eroding.
    slice_parapet = (slice_outer.*slice_Z_all) > (height_roof+para_minheight);
    slice_parapet = imdilate(slice_parapet,strel('disk',6));
    slice_parapet = imerode(slice_parapet,strel('disk',5));   % `slice_parapet` is somewhat larger than before this way.
    slice_parapet = slice_parapet & slice_all;   % Makes sure that slice_parapet does not reach beyond roof.

    % Get final parapet height!
    if sum(slice_parapet,'all')>0
        height_parapet = nansum(slice_parapet.*slice_Z_all,'all') / sum(slice_parapet,'all') - height_roof;
    else
        height_parapet = NaN;
    end

    % CLEANUP: Delete any blob from final no-parapet array that's way too tiny.
    slice_all_noparapet = bwareaopen(slice_all.*(~slice_parapet),ceil(sum(slice_all,'all')*0.05),8);
end