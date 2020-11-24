function names = file_exists(folder,extension)
%FILE_EXISTS returns a string array of filenames (EXCLUDING extensions)
%that contain extension `extension` in folder `folder`.
%
%   Parameters:
%   -----------
%   folder : string
%       Folder name (e.g. "subfolder", "subfolder\", or "").
%   extension : string
%       File extension (e.g. ".png").
%
%   Outputs:
%   --------
%   names : string
%       Array of filenames ending with `extension` in `folder`.

    file_list = dir(fullfile(folder, "*"));
    folder_contents = {file_list.name};
    filenames = string(folder_contents');
    file_exists_idx = contains(filenames,extension);

    names = filenames(file_exists_idx);
    names = erase(names,extension);
end