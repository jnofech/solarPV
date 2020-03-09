function isroof = find_roofs(issurface, roof_mindiameter)
%FIND_SURFACES() identifies roofs (i.e. any surface large enough to fit a
%   tiny solar panel under California regulations regarding edge distance).
%
%   Parameters:
%   -----------
%   issurface : logical (2D array)
%       Indicates any structure that is not the ground, an edge, or tiny.
%   roof_mindiameter : double
%       Minimum diameter of the shortest side of the rooftop, in pixels.
%       (This should be related to the shortest side of the solar PV.)
%
%   Returns:
%   --------
%   isroof : logical (2D array)
%       Indicates any surface large enough to hold a tiny solar panel.
    
    % Create `Z_isroof`-- Identify blobs that can fit a solar panel
    blobs = blob_clean(issurface,roof_mindiameter);
    isroof = logical(blobs);
end

