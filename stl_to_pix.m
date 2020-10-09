function varargout = stl_to_pix(fname,length,theta_CCW)
%STL_TO_PIX Reads STL with `stlread.m`, converting units to pixels.
%   Parameters:
%   -----------
%   fname : string
%       Name of STL file (with extension), including folder.
%   length : integer
%       Desired pixel length of longer side of model. Higher numbers mean
%       better resolution, but longer computing times.
%
%   Outputs:
%   --------
%   `p = stl_to_pix(fname,length)`:
%       `p` is a structure containing faces, vertices, and normals of mesh.
%   `[f,v] = stl_to_pix(fname,length)`:
%       `f`,`v` are faces and vertices of mesh.
%   `[f,v,n] = stl_to_pix(fname,length)`:
%       `f`,`v`,`n` are faces, vertices, and normals of mesh.
    % 3D rotation matrix
    theta_CCW = -theta_CCW;     % Swap direction (X and Y are swapped in matlab?)
    Rz = [  cosd(theta_CCW), -sind(theta_CCW), 0;
            sind(theta_CCW),  cosd(theta_CCW), 0;
                          0,                0, 1];
    p = stlread(fname);
    % Flip y-axis, to match Matlab's imshow/imagesc defaults. Flip normals
    % to match.
    p.vertices(:,2) = -1*p.vertices(:,2);
    p.normals(:,2) = -1*p.normals(:,2);
    % Shift all vertices such that none are negative
    p.vertices(:,1) = p.vertices(:,1) - min(p.vertices(:,1));
    p.vertices(:,2) = p.vertices(:,2) - min(p.vertices(:,2));
    p.vertices(:,3) = p.vertices(:,3) - min(p.vertices(:,3));
    % Scale to desired length
    p.vertices = p.vertices / max(p.vertices(:,1:2),[],'all') * floor(length-1); % Only consider x,y lengths. Not height!
    p.vertices = p.vertices+1; % One last nudge, so that no vertices are on x,y,z=0 (which cannot be indexed)
    
    % Rotate!
    p.vertices = (Rz * p.vertices')';
    p.normals = (Rz * p.normals')';
    % Shift all vertices once more, such that none are negative
    p.vertices(:,1) = p.vertices(:,1) - min(p.vertices(:,1));
    p.vertices(:,2) = p.vertices(:,2) - min(p.vertices(:,2));
    p.vertices(:,3) = p.vertices(:,3) - min(p.vertices(:,3));
    p.vertices = p.vertices+1; % One last nudge, so that no vertices are on x,y,z=0 (which cannot be indexed)

%     varargout = cell(1,nargout);
    switch nargout        
        case 2
            varargout{1} = p.faces;
            varargout{2} = p.vertices;
        case 3
            varargout{1} = p.faces;
            varargout{2} = p.vertices;
            varargout{3} = p.normals;
        otherwise
            varargout{1} = struct('faces',p.faces,'vertices',p.vertices,'normals',p.normals);
    end
end

