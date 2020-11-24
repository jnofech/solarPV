function automate_model_generation(building_URLs,varargin)
% AUTOMATE_MODEL_GENERATION generates 3D models at the user-specified
% location(s) using RenderDoc, Chrome/Chromium, and Blender.
% WARNING: You will NOT be able to work in the background while this is 
% running, as it makes use of direct mouse/keyboard inputs.
%
% Parameters:
% -----------
% building_URLs : cell array
%   An Nx2 cell array containing URLs and names, where N is the number of
%   buildings. The URL refers to a direct link to each building's position
%   (with zoom) on Google Maps. The `name` is simply the name of the
%   building, preferably without spaces or special characters.
%   Example:
%       building_URLs = {"https://www.google.com/maps/@?api=1&map_action=map&center=53.5281333,-113.5256635&zoom=19&basemap=satellite", "TEST";
%                        "https://www.Google.com/maps/@35.3626661,138.7283781,3000m/data=!3m1!1e3",                                     "MtFuji"};
%
% `path_meshes` ("<main path>+\meshes\") : string
%   FULL path to output 3D meshes ("<name>.stl").
%   Blender project files are also saved in this path ("<name>.blend").
%
% `path_images` ("<main path>+\images\") : string
%   FULL path to output 2D satellite photos from Google Maps ("<name>.png").
%   This function automatically saves in .png format, but .jpg is also
%   readable by the main rooftop segmentation script.
%
% `path_RDC` ("<main path>+\meshes\") : string
%   FULL path to savefile location for RenderDoc captures ("<name>.rdc").
%
% `path_buttons` ("\output\temp_buttons\") : string
%   Path to savefile location for certain button locations. These are saved
%   after being generated for the first time.
%
% `timer_mult` (1) : double
%   Multiplier for the timed "pauses" that this function takes to allow
%   things to load. Larger values lead to safer and more consistent
%   performance, but of course take longer.
%
%
% Outputs:
% For each input `building_URLs` of {"<URL>","<name>"}, this code produces:
% --------
% "<path_RDC>\<name>.rdc"       : RenderDoc GPU capture
% "<path_images>\<name>.png"    : 2D satellite photo
% "<path_images>\<name>.mat"    : 2D satellite photo height in metres
% "<path_meshes>\<name>.stl"    : 3D model
% "<path_meshes>\<name>.blend"  : Blender project

    % Set default values
    global timer_mult
    timer_mult = 1;
%     path_buttons = "output/temp_buttons/";
    path_buttons = pwd+"\output\";
    % path_RDC     = "C:\Users\Win10\Desktop\Building Captures";
%     path_RDC     = pwd+"\RDC\";
    path_RDC     = pwd+"\meshes\";
    path_meshes  = pwd+"\meshes\";
    path_images  = pwd+"\images\";

    if nargin >= 2
        for idx = 1:2:length(varargin)
            switch lower(varargin{idx})
                case "path_buttons"
                    path_buttons = varargin{idx+1};
                case "path_rdc"
                    path_RDC = varargin{idx+1};
                case "path_meshes"
                    path_meshes = varargin{idx+1};
                case "path_images"
                    path_images = varargin{idx+1};
                case "timer_mult"
                    timer_mult = varargin{idx+1};
            end
        end
    end
    % Convert all paths to Strings
    path_buttons = string(path_buttons);
    path_RDC = string(path_RDC);
    path_meshes = string(path_meshes);
    path_images = string(path_images);
    
    n_buildings = size(building_URLs,1);

    % FIRST: Check if user wants to overwrite meshes+images
    names_stl  = file_exists(path_meshes,".stl");
    names_png  = file_exists(path_images,".png");
    names_jpg  = file_exists(path_images,".jpg");
    names_jpeg = file_exists(path_images,".jpeg");
    names_bmp  = file_exists(path_images,".bmp");
    idx_skip = [];  % Index of buildings that are SKIPPED because user did not want to overwrite.

    for i = 1:size(building_URLs,1)
        name = building_URLs{i,2};
        if ismember(name,names_stl) || (ismember(name,names_png) || ismember(name,names_jpg) || ismember(name,names_jpeg) || ismember(name,names_bmp))
            % Decide whether to overwrite.
            output = questdlg("3D mesh and/or satellite photo for "+name+" already exist. Overwrite?","Confirm overwrite","Yes","No","Yes");
            if isequal(output,'Yes')
                % Do nothing
            else
                idx_skip = [idx_skip,i];
            end
        end
    end

    % Set file names for button locations
    file_RenderDoc_pid      = "button_RenderDoc_pid.mat";
    file_RenderDoc_timer    = "button_RenderDoc_timer.mat";
    file_Chrome_3D_switcher = "button_Chrome_3D_switcher.mat";
    file_Blender_path_bar   = "button_Blender_path_bar.mat";
    file_Blender_filename_bar = "button_Blender_filename_bar.mat";
    file_Chrome_map_window  = "button_Chrome_map_window.mat";

    % Begin script
%     clc
    close all
    output = questdlg("Please open Blender and close any instances of Google Chrome before continuing.","Confirm: open Blender","Continue","Cancel","Continue");
    if isequal(output,'Cancel')
        error("Canceled.");
    end


    % Initialize external apps
    if cmd_isopen("chromium") && cmd_isopen("renderdoc")
        rd_injected = 1;
    elseif ~(cmd_isopen("chromium") && cmd_isopen("renderdoc"))
        % Open Chromium & RenderDoc anew!
        rd_injected = 0;
        disp("Initializing: RenderDoc, Chromium");
        cmd_open("RenderDoc",0.5);
        cmd_open("Chromium",0.5);

        % Wait for external apps to open.
        % cmd_wait_until_open("blender",0);
        % disp("... Blender is open.")
        cmd_wait_until_open("renderdoc",0);
        disp("... RD is open.")
        cmd_wait_until_open("chromium",0);
        disp("... Chrome/chromium is open.")

        % Let user determine when everything's ready to begin.
        output = questdlg("Click 'Continue' when RenderDoc, Blender, and Chromium are ready to use. Be sure to close RenderDoc splash window first.","Confirm: Ready?","Continue","Cancel","Continue");
        if isequal(output,'Continue')
            % Do nothing
            pause(timer_mult*1);
        else
            error("Canceled.");
        end
    else
        error("Close RenderDoc and Chromium and try again.");
    end

    % Initialize Java Robot
    % (!) Be aware that I should import every time I'm about to make an input!
    import java.awt.Robot;
    mouse = Robot;
    % import java.awt.event.*;
    import java.awt.event.KeyEvent;
    import java.awt.event.InputEvent;



    if ~rd_injected
        % Inject Chromium into RenderDoc!
        % ~~ CHROMIUM ~~
    %         chrom_mode = "pid";
            chrom_mode = "chromium gpu";
            if isequal(chrom_mode,"pid")
                % Get user to copy Chromium GPU 'pid'
                cmd_switchto("chromium",0.5);
                pause(timer_mult*0.2);
                pid = inputdlg("From 'Chromium GPU' popup, input GPU 'pid':","'pid' confirmation" );
                pid = str2double(pid{1});
                while isnan(pid)
                    pid = inputdlg("Invalid 'pid'. From 'Chromium GPU' popup, input GPU 'pid':","'pid' confirmation" );
                    pid = str2double(pid{1});
                end
                pid = num2str(pid); % Can loop through each character!
                warning("There may be multiple processes containing PID of "+string(pid)+"!");
            elseif isequal(chrom_mode,"chromium gpu")
                pid = 'CHROMIUM GPU';
            else
                error("Invalid `chrom_mode`.");
            end

        % ~~ RENDERDOC ~~
            % RenderDoc: Put into injection mode.
            cmd_switchto("renderdoc",0);

            % Tap: Win+Up (Maximize)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_UP)
            mouse.keyRelease(KeyEvent.VK_WINDOWS)
            pause(timer_mult*0.1)

            % Tap: Alt+(F+I) + Enter
            mouse.keyPress(KeyEvent.VK_ALT)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_F)
                mouse.keyRelease(KeyEvent.VK_F)
                mouse.keyPress(KeyEvent.VK_I)
                mouse.keyRelease(KeyEvent.VK_I)
            mouse.keyRelease(KeyEvent.VK_ALT)
            pause(timer_mult*0.1)
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*1)

            % Get screen size
            screenSize = get(0, 'screensize');
            left = screenSize(1);
            top  = screenSize(2);
            width  = screenSize(3);
            height = screenSize(4);

            % Get button (i.e. text box) position!
            if isfile(path_buttons+file_RenderDoc_pid)
                load(path_buttons+file_RenderDoc_pid,"filter_x","filter_y");
            else
                % User inputs button position for the first time.
                msg = "Below, click on position of 'Filter process list by PID or name' text box (shown above).";
                snippet = imread('RenderDoc_textbox.png');
                [filter_x,filter_y] = ui_getButton(msg,snippet,10);
                save(path_buttons+file_RenderDoc_pid,"filter_x","filter_y");
                pause(1);
            end

            % Click on dialog box location, and type in `pid`!
            move(mouse,[filter_x,filter_y],'mode','absolute','mult',2);
            tap(mouse,1);
            type_keys(mouse,pid);
            % Click above dialog box
            pause(timer_mult*0.02);
            move(mouse,[0,-50],'mode','relative','mult',2);
            pause(timer_mult*0.5);
            tap(mouse,1);
            % Select process (only 1 should be available!), and activate it
            import java.awt.event.KeyEvent;
            for ii = 1:50
                mouse.keyPress(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_UP)
                pause(timer_mult*0.001)
            end
            pause(timer_mult*0.5);
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*5);

        % ~~ CHROMIUM ~~
            % Chromium: Close PID dialog.
            cmd_switchto("chromium",0);
            pause(timer_mult*0.5);
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*2);

        rd_injected = true;
    else
        % Do nothing
    end

    % Loop through buildings!
    for i = 1:n_buildings
%     for i = 1
        URL = building_URLs{i,1};
        building_name = string(building_URLs{i,2});

        % Skip buildings if you didn't want to overwrite!
        if ismember(i,idx_skip)
            continue
        end

            % Get screen size
            screenSize = get(0, 'screensize');
            left = screenSize(1);
            top  = screenSize(2);
            width  = screenSize(3);
            height = screenSize(4);

            % Chromium: Go to browser. Maximize, and begin entering URL.
            cmd_switchto("chromium",0.3);
            % Tap: Win+Up (Maximize)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_UP)
            mouse.keyRelease(KeyEvent.VK_WINDOWS)
            pause(timer_mult*0.3)
            % Tap: CTRL+(L)
            mouse.keyPress(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_L)
                mouse.keyRelease(KeyEvent.VK_L)
            mouse.keyRelease(KeyEvent.VK_CONTROL)
            pause(timer_mult*0.1)
            % Type address!
            type_keys(mouse,URL);
            % Tap: ENTER
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*0.3)

            % Wait for GMaps to load.
            warning("Waiting times for Google Maps loading can vary depending on internet speed. Consider checking bandwidth usage?");
            pause(timer_mult*2);

            % Find button location for 3D switcher
            if isfile(path_buttons+file_Chrome_3D_switcher)
                load(path_buttons+file_Chrome_3D_switcher,"switch_x","switch_y");
            else
                % User inputs button position for the first time.
                msg = "Below, click on position of 'Enable/Disable Globe View' button (shown above).";
                snippet = imread('Google_Set3D.png');
                [switch_x,switch_y] = ui_getButton(msg,snippet,6);
                save(path_buttons+file_Chrome_3D_switcher,"switch_x","switch_y");
                pause(1);
            end

            % Ask user whether it's currently in 2D or 3D mode
    %         output = questdlg("Is Google Maps displaying in 2D or 3D?","Google Maps 2D/3D","2D","3D","2D");
            output = "2D"; warning("CHROMIUM: Assuming that GMaps begins in 2D mode.");

            % Get height from address bar. If this doesn't work, retry
            % a few times until it does. Raise error if it fails too many
            % times.
            URL_parse_failed = false;
            URL_parse_max_duration = 20*timer_mult;
            GMaps_zoom_applied = false;     % Tracks whether zoom was modified yet
            tic
            while toc < URL_parse_max_duration
                % Copy address
                % Tap: CTRL+(L+C)
                pause(timer_mult*0.0)
                mouse.keyPress(KeyEvent.VK_CONTROL)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_L)
                    mouse.keyRelease(KeyEvent.VK_L)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_C)
                    mouse.keyRelease(KeyEvent.VK_C)
                mouse.keyRelease(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)

                % Get latitude from address
                expression = '@[-1234567890.]+,';
                url_active = clipboard('paste');
                lat = regexp(url_active,expression,'match');
                try
                    % Read latitude from URL
                    lat = erase(lat{:},'@');
                    lat = erase(lat,',');
                    url_active_alt = regexprep(url_active,expression,'_LATITUDE_,');
                    % Read longitude from modified URL
                    expression = '_LATITUDE_,[-1234567890.]+,';
                    lon = regexp(url_active_alt,expression,'match');
                    lon = erase(lon{:},'_LATITUDE_,');
                    lon = erase(lon,',');
                    url_active_alt = regexprep(url_active_alt,expression,'_LONGITUDE_,');
                    % Read height (m) from modified URL
                    expression = '_LONGITUDE_,[1234567890.]+m';
                    cmap_height = regexp(url_active_alt,expression,'match');
                    cmap_height = erase(cmap_height{:},'_LONGITUDE_,');
                    cmap_height = erase(cmap_height,'m');
                    save(path_images+building_name+".mat","cmap_height");
                    URL_parse_failed = false;
                    
                    % If successful so far: Move to center of screen, and
                    % zoom in once!
                    if GMaps_zoom_applied == false
                        % Re-center mouse (to keep building centered)
                        move(mouse,[width/2,height/2],"absolute");   
                        % Zoom out once
%                         mouse.mouseWheel(-1);
                        % (REMOVED since zooming too far in causes terrain
                        % to stop rendering properly.)
                        pause(1);
                        GMaps_zoom_applied = true;
                    else
                        break
                    end
                catch
                    if toc > URL_parse_max_duration/2
                        disp("Failed to parse latitude, longitude, and/or height from URL. Retrying.");
                    end
                    URL_parse_failed = true;
                end
            end
            if URL_parse_failed
                error("Unable to parse latitude, longitude, and/or height from URL '"+url_active+"'."+newline ...
                      +"Check that URL is in format 'https://www.google.com/maps/@<LATITUDE>,<LONGITUDE>,<HEIGHT>m/data=<...>'. If not, then `automate_model_generation.m` will have to be updated accordingly.");
            end
            
            url_active_alt = regexprep(url_active_alt,expression,'_LATITUDE_,');

            % Take fullscreen snapshot in 2D
            img = takeScreenshot(left,top,width,height);

            % Ask user to click on upper+lower boundaries of GMaps window (saved permanently; used to get pix2m for each unique building)
            if isfile(path_buttons+file_Chrome_map_window)
                load(path_buttons+file_Chrome_map_window,"corners_x","corners_y");
            else
                [corners_x,corners_y] = getBorders(left,top,width,height, "Select corners of full Google Maps interface. Be as accurate as possible!");
                save(path_buttons+file_Chrome_map_window,"corners_x","corners_y");
                pause(1);
            end
            ydist_pix = abs(corners_y(1)-corners_y(2));     % Vertical pixel distance of GMaps UI
            % Save satellite photo for later use!
            img_cropped = imcrop(img,[corners_x(1), corners_y(1), corners_x(2)-corners_x(1), corners_y(2)-corners_y(1)]);
            imwrite(img_cropped,path_images+building_name+".png"); close;

            % Toggle to 3D
            if isequal(output,"2D")
                cmd_switchto("chromium",1);
                move(mouse,[switch_x,switch_y],'mode','absolute','mult',2);
                tap(mouse,1);
                % Give a bit of time to load
                pause(timer_mult*3);
            else
                % Do nothing
            end

        % ~~ RENDERDOC ~~
            % Switch to RenderDoc
            cmd_switchto("RenderDoc",0.2);
            % Tap: Win+Up (Maximize)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_UP)
            mouse.keyRelease(KeyEvent.VK_WINDOWS)
            pause(timer_mult*0.7)

            % Ask user for "0 secs" button location
            if isfile(path_buttons+file_RenderDoc_timer)
                load(path_buttons+file_RenderDoc_timer,"timer_x","timer_y");
            else
                % User inputs button position for the first time.
                pause(timer_mult*1);
                msg = "Below, click on position of '0 secs' box (shown above).";
                snippet = imread('RenderDoc_timer.png');
                [timer_x,timer_y] = ui_getButton(msg,snippet,6);
                save(path_buttons+file_RenderDoc_timer,"timer_x","timer_y");
                pause(1);
            end

            % Triple-click on "0 secs" button (backspace doesn't work)
            move(mouse,[timer_x,timer_y],'mode','absolute','mult',2);
            tripletap(mouse,1);
            pause(timer_mult*0.2);

            % Set timer (4 secs default)
            import java.awt.event.KeyEvent;
            type_keys(mouse,num2str(ceil(4*timer_mult)));

            % Move a few pixels to the right; click
            move(mouse,[timer_x,timer_y],'mode','absolute','mult',2);
            move(mouse,[75,0],'mode','relative','mult',2);
            tap(mouse,1);
            move(mouse,[-75,0],'mode','relative','mult',1);
            pause(timer_mult*0.1);

            % Switch to Chromium
            cmd_switchto("chromium",1);

            % Move to middle of screen; drag around in r=10pix circles for 4 seconds
            circ_radius = 30;   % in pixels
            circ_rps = 3;       % Rotations per second. Becomes inaccurate above 3
            if cmap_height < 120
                % Height is too low, and the ground won't be moving much
                % relative to the buildings due to GMaps' built-in 
                % perspective distortion.
                % Compensate for this with wider mouse rotations.
                circ_radius = circ_radius*1;
            end
            move(mouse,[width/2+circ_radius,height/2],'mode','absolute','mult',1);
            press(mouse,1);
            tic
            while toc < ceil(4*timer_mult)
                % Drag cursor in a circle of fluctuating radius.
                % (Parts of the ground may be missing if the mouse moves
                % in too consistent of a pattern.)
                rad = circ_radius + (circ_radius/3)*sin(2*pi*toc*2);
                x_offset = rad*cos(2*pi*toc*circ_rps);
                y_offset = rad*sin(2*pi*toc*circ_rps);
                teleport(mouse,[width/2+x_offset,height/2+y_offset],'mode','absolute','mult',1);
            end
            pause(timer_mult*0.1);
            release(mouse,1);

            % Switch back to 2D
            cmd_switchto("chromium",1);
            move(mouse,[switch_x,switch_y],'mode','absolute','mult',2);
            tap(mouse,1);
            pause(timer_mult*1);

            % Switch to RenderDoc
            cmd_switchto("renderdoc",1);

            % Move down-left of timer button; tap M2 + S
            % (SAVE capture!)
            move(mouse,[timer_x,timer_y],'mode','absolute','mult',2);
            move(mouse,[-150,150],'mode','relative','mult',2);
            pause(timer_mult*0.1);
            tap(mouse,3);
            pause(timer_mult*0.1);
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_S)
            mouse.keyRelease(KeyEvent.VK_S)
            pause(timer_mult*1);

            % Press "Tab" 6 times, then ENTER
            % (select path)
            tabs_to_path = 6;
            import java.awt.event.KeyEvent;
            for ii = 1:tabs_to_path
                mouse.keyPress(KeyEvent.VK_TAB)
                mouse.keyRelease(KeyEvent.VK_TAB)
                pause(timer_mult*0.02);
            end
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*0.1)

            % Type in RDC path, then tap ENTER
    %         path_RDC = "C:\Users\Win10\Desktop\Building Captures";
            type_keys(mouse,path_RDC);
            pause(timer_mult*0.1)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*0.3)

            % Press "Shift+TAB" 8 times
            % (The extra 2 times are because of new buttons that activated
            % next to path box.)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_SHIFT)
                for ii = 1:(tabs_to_path+2)
                    mouse.keyPress(KeyEvent.VK_TAB)
                    mouse.keyRelease(KeyEvent.VK_TAB)
                    pause(timer_mult*0.02);
                end
            mouse.keyRelease(KeyEvent.VK_SHIFT)
            pause(timer_mult*0.3)

            % Type in building name, then save
    %         building_name = "TEST";
    %         building_name = "MtFuji";
            type_keys(mouse,building_name);
            overwrite = true;
            pause(timer_mult*0.5);
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*0.2);
            if overwrite
                % Tap "yes" at dialog box. (No effect on files that don't
                % previously exist.)
                mouse.keyPress(KeyEvent.VK_Y)
                mouse.keyRelease(KeyEvent.VK_Y)
            else
                % Tap "no" at dialog box. (No effect on files that don't
                % previously exist.)
                mouse.keyPress(KeyEvent.VK_N)
                mouse.keyRelease(KeyEvent.VK_N)
            end
            pause(timer_mult*1);

            % Move down-left of timer button; tap M2 + D + Y
            % (DELETE capture!)
            move(mouse,[timer_x,timer_y],'mode','absolute','mult',2);
            move(mouse,[-150,150],'mode','relative','mult',2);
            pause(timer_mult*0.1);
            tap(mouse,3);
            pause(timer_mult*0.1);
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_D)
            mouse.keyRelease(KeyEvent.VK_D)
            pause(timer_mult*0.1);
            mouse.keyPress(KeyEvent.VK_Y)
            mouse.keyRelease(KeyEvent.VK_Y)
            pause(timer_mult*1);

    %     % ~~ BLENDER ~~
    %         cmd_open("blender",1);    % DOES NOT WORK! Blender needs to be opened beforehand.
    %         pause(timer_mult*5);
            cmd_switchto("blender",1);
            % Tap: Win+Down, Win+Up (Maximize)
            import java.awt.event.KeyEvent;
            mouse.keyPress(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_UP)
            mouse.keyRelease(KeyEvent.VK_WINDOWS)
            pause(timer_mult*0.3)
            % Tap: Esc (close splash window, if any)
            mouse.keyPress(KeyEvent.VK_ESCAPE)
            mouse.keyRelease(KeyEvent.VK_ESCAPE)
            pause(timer_mult*0.5);
            % Re-center mouse (important in Blender)
            move(mouse,[width/2,height/2],"absolute");      
            % Delete everything! (A, X, Enter)
            % ("VK_DELETE" does not work properly in Blender for some reason)
            mouse.keyPress(KeyEvent.VK_A)
            mouse.keyRelease(KeyEvent.VK_A)
            pause(timer_mult*0.3);
            mouse.keyPress(KeyEvent.VK_X)
            mouse.keyRelease(KeyEvent.VK_X)
            pause(timer_mult*0.1);
            mouse.keyPress(KeyEvent.VK_ENTER)
            mouse.keyRelease(KeyEvent.VK_ENTER)
            pause(timer_mult*0.3);

            % Tap: F4, I, R (import .RDC file)
            mouse.keyPress(KeyEvent.VK_F4)
            mouse.keyRelease(KeyEvent.VK_F4)
            pause(timer_mult*0.1)
            mouse.keyPress(KeyEvent.VK_I)
            mouse.keyRelease(KeyEvent.VK_I)
            pause(timer_mult*0.1)
            mouse.keyPress(KeyEvent.VK_R)
            mouse.keyRelease(KeyEvent.VK_R)
            pause(timer_mult*0.5)

                % Tap: Win+Up (Maximize)
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_WINDOWS)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_UP)
                    mouse.keyRelease(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.3)

                % Ask user for "Path" input location
                if isfile(path_buttons+file_Blender_path_bar)
                    load(path_buttons+file_Blender_path_bar,"path_x","path_y");
                else
                    % User inputs button position for the first time.
                    msg = "Below, click on position of 'path' text box (at top).";
                    snippet = [];
                    [path_x,path_y] = ui_getButton(msg,snippet,6);
                    save(path_buttons+file_Blender_path_bar,"path_x","path_y");
                    pause(1);
                end
                % Ask user for "Filename" input location
                if isfile(path_buttons+file_Blender_filename_bar)
                    load(path_buttons+file_Blender_filename_bar,"filename_x","filename_y");
                else
                    % User inputs button position for the first time.
                    msg = "Below, click on position of 'filename' text box (at bottom).";
                    snippet = [];
                    [filename_x,filename_y] = ui_getButton(msg,snippet,6);
                    save(path_buttons+file_Blender_filename_bar,"filename_x","filename_y");
                    pause(1);
                end
                % Move to path textbox location, and click on it
                move(mouse,[path_x,path_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);        

                % CTRL+A, type path, Enter
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_CONTROL)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_A)
                    mouse.keyRelease(KeyEvent.VK_A)
                mouse.keyRelease(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)
                type_keys(mouse,path_RDC);
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)

                % Move to path filename textbox location, and click on it
                move(mouse,[filename_x,filename_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);  
                type_keys(mouse,building_name+".rdc");
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                % Re-center mouse (important in Blender)
                move(mouse,[width/2,height/2],"absolute");

            % Wait for Blender to fully load the .STL
            cmd_wait_for_response("blender",2);

            % Tap: F4, E, T (export .STL file)
            mouse.keyPress(KeyEvent.VK_F4)
            mouse.keyRelease(KeyEvent.VK_F4)
            pause(timer_mult*0.1)
            mouse.keyPress(KeyEvent.VK_E)
            mouse.keyRelease(KeyEvent.VK_E)
            pause(timer_mult*0.1)
            mouse.keyPress(KeyEvent.VK_T)
            mouse.keyRelease(KeyEvent.VK_T)
            pause(timer_mult*0.5)

                % Tap: Win+Up (Maximize)
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_WINDOWS)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_UP)
                    mouse.keyRelease(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.3)

                % Move to path textbox location, and click on it
                move(mouse,[path_x,path_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);        

                % CTRL+A, type path, Enter
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_CONTROL)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_A)
                    mouse.keyRelease(KeyEvent.VK_A)
                mouse.keyRelease(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)
                type_keys(mouse,path_meshes);
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.1);

                % Move to path filename textbox location, and click on it
                move(mouse,[filename_x,filename_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);  
    %             type_keys(mouse,building_name+"_entire.stl");
                type_keys(mouse,building_name+".stl");
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                % Re-center mouse (important in Blender)
                move(mouse,[width/2,height/2],"absolute");
                tap(mouse,1);
                pause(timer_mult*0.5)


            % Tap: CTRL+SHIFT+S (save as .blend file)
            mouse.keyPress(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)
                mouse.keyPress(KeyEvent.VK_SHIFT)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_S)
                    mouse.keyRelease(KeyEvent.VK_S)
                    pause(timer_mult*0.1)
                mouse.keyRelease(KeyEvent.VK_SHIFT)
                pause(timer_mult*0.1)
            mouse.keyRelease(KeyEvent.VK_CONTROL)
            pause(timer_mult*0.5)

                % Tap: Win+Up (Maximize)
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_WINDOWS)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_UP)
                    mouse.keyRelease(KeyEvent.VK_UP)
                mouse.keyRelease(KeyEvent.VK_WINDOWS)
                pause(timer_mult*0.3)

                % Move to path textbox location, and click on it
                move(mouse,[path_x,path_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);        

                % CTRL+A, type path, Enter
                import java.awt.event.KeyEvent;
                mouse.keyPress(KeyEvent.VK_CONTROL)
                    pause(timer_mult*0.1)
                    mouse.keyPress(KeyEvent.VK_A)
                    mouse.keyRelease(KeyEvent.VK_A)
                mouse.keyRelease(KeyEvent.VK_CONTROL)
                pause(timer_mult*0.1)
                type_keys(mouse,path_meshes);
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)

                % Move to path filename textbox location, and click on it
                move(mouse,[filename_x,filename_y],'mode','absolute','mult',2);
                pause(timer_mult*0.1);
                tap(mouse,1);
                pause(timer_mult*0.1);  
                type_keys(mouse,building_name+".blend");
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                mouse.keyPress(KeyEvent.VK_ENTER)
                mouse.keyRelease(KeyEvent.VK_ENTER)
                pause(timer_mult*0.3)
                % Re-center mouse (important in Blender)
                move(mouse,[width/2,height/2],"absolute");

            % Tap: F4, N, G, D (create new project, so window name is "blender" again)
            mouse.keyPress(KeyEvent.VK_F4)
            mouse.keyRelease(KeyEvent.VK_F4)
            pause(timer_mult*0.1);
            mouse.keyPress(KeyEvent.VK_N)
            mouse.keyRelease(KeyEvent.VK_N)
            pause(timer_mult*0.1);
            mouse.keyPress(KeyEvent.VK_G)
            mouse.keyRelease(KeyEvent.VK_G)
            pause(timer_mult*0.1);
            mouse.keyPress(KeyEvent.VK_D)
            mouse.keyRelease(KeyEvent.VK_D)
            pause(timer_mult*0.3);

    end
end

% Mouse Functions

function coords = locate()
% Returns locations of cursor's current position, in coordinates used by
% Java robot. (i.e. top left is [1,1], and bottom right is
% [screenwidth,screenheight].)
    screenSize  = get(0, 'screensize');
    coords      = get(0, 'PointerLocation');
    coords(2) = screenSize(4) - coords(2) + 1;
end

function tap(mouse,button)
% Taps specified mouse button.
    global timer_mult
    press(mouse,button);
    pause(timer_mult*0.001);
    release(mouse,button);
end

function doubletap(mouse,button)
% Double-taps specified mouse button.
    global timer_mult
    tap(mouse,button);
    pause(timer_mult*0.01);
    tap(mouse,button);
end

function tripletap(mouse,button)
% Triple-taps specified mouse button.
    global timer_mult
    tap(mouse,button);
    pause(timer_mult*0.01);
    tap(mouse,button);
    pause(timer_mult*0.01);
    tap(mouse,button);
end

function drag(mouse,button,coords,varargin)
% Clicks+drags mouse button from current location to [x,y] `coords`.
% The `mode` ("rel" or "abs") determines whether the coordinates are
% relative to cursor location, or absolute.
    global timer_mult
    if nargin-3 == 0
        mode = "relative";
    elseif nargin-3 == 1
        mode = varargin(1);
        mode = mode{:};
    else
        error("Too many input arguments.");
    end
    press(mouse,button)
    move(mouse,coords,mode)
    release(mouse,button)
end

function move(mouse,coords,varargin)
% Moves mouse button from current location to specified [x,y] `coords`.
% The `mode` ("rel" or "abs") determines whether the coordinates are
% relative to cursor location, or absolute.
    % Read inputs
    switch length(varargin)
        case 0
            mode = "relative";
            step_multiplier = 1;
        case 1
            mode = varargin{1};
            step_multiplier = 1;
        otherwise
            mode = "relative";
            step_multiplier = 1;
            for i = 1:length(varargin)
                if isequal(lower(varargin{i}),'mode')
                    mode = varargin{i+1};
                elseif isequal(lower(varargin{i}),'step_multiplier') || ...
                       isequal(lower(varargin{i}),'step_mult') || ...
                       isequal(lower(varargin{i}),'time_mult') || ...
                       isequal(lower(varargin{i}),'time') || ...
                       isequal(lower(varargin{i}),'timemult') || ...
                       isequal(lower(varargin{i}),'mult') || ...
                       isequal(lower(varargin{i}),'slow')
                    step_multiplier = varargin{i+1};
                end
            end
    end
    
    % Determine if specified coordinates are relative or absolute.
    location = locate();
    if isequal(mode,1) || isequal(mode,"relative") || isequal(mode,"rel")
        % `coords` is relative to cursor location.
        truecoords = coords + location;
    elseif isequal(mode,2) || isequal(mode,"absolute") || isequal(mode,"abso")
        % `coords` is absolute cursor location.
        truecoords = coords;
    else
        error("Invalid mouse dragging mode.");
    end
    
    % Move mouse!
    time_per_step = 0.0005;  % In seconds.
    steps = max(abs(truecoords - location)) / 5 * step_multiplier;
    dx = (truecoords(1) - location(1)) / steps;
    dy = (truecoords(2) - location(2)) / steps;
    for i = 1:steps
        mouse.mouseMove(location(1) + dx*i, location(2) + dy*i);
        pause(time_per_step);
    end
    mouse.mouseMove(truecoords(1),truecoords(2));
end

function teleport(mouse,coords,varargin)
% Moves mouse button from current location to specified [x,y] `coords`.
% The `mode` ("rel" or "abs") determines whether the coordinates are
% relative to cursor location, or absolute.
% Moves in a single step.
    % Read inputs
    switch length(varargin)
        case 0
            mode = "relative";
            step_multiplier = 1;
        case 1
            mode = varargin{1};
            step_multiplier = 1;
        otherwise
            mode = "relative";
            step_multiplier = 1;
            for i = 1:length(varargin)
                if isequal(lower(varargin{i}),'mode')
                    mode = varargin{i+1};
                elseif isequal(lower(varargin{i}),'step_multiplier') || ...
                       isequal(lower(varargin{i}),'step_mult') || ...
                       isequal(lower(varargin{i}),'time_mult') || ...
                       isequal(lower(varargin{i}),'time') || ...
                       isequal(lower(varargin{i}),'timemult') || ...
                       isequal(lower(varargin{i}),'mult') || ...
                       isequal(lower(varargin{i}),'slow')
                    step_multiplier = varargin{i+1};
                end
            end
    end
    
    % Determine if specified coordinates are relative or absolute.
    location = locate();
    if isequal(mode,1) || isequal(mode,"relative") || isequal(mode,"rel")
        % `coords` is relative to cursor location.
        truecoords = coords + location;
    elseif isequal(mode,2) || isequal(mode,"absolute") || isequal(mode,"abso")
        % `coords` is absolute cursor location.
        truecoords = coords;
    else
        error("Invalid mouse dragging mode.");
    end
    
    % Move mouse!
    time_per_step = 0.00000;  % In seconds.
%     steps = max(abs(truecoords - location)) / 0.5 * step_multiplier;
    steps = 1;
    dx = (truecoords(1) - location(1)) / steps;
    dy = (truecoords(2) - location(2)) / steps;
    for i = 1:steps
        mouse.mouseMove(location(1) + dx*i, location(2) + dy*i);
        pause(time_per_step);
    end
    mouse.mouseMove(truecoords(1),truecoords(2));
end

function press(mouse,button)
% Presses specified mouse button.
    import java.awt.Robot;
    % import java.awt.event.*;
    import java.awt.event.KeyEvent;
    import java.awt.event.InputEvent;
    switch button
        case 1
            mouse.mousePress(InputEvent.BUTTON1_MASK);
        case 2
            mouse.mousePress(InputEvent.BUTTON2_MASK);
        case 3
            mouse.mousePress(InputEvent.BUTTON3_MASK);
        case 4
            mouse.mousePress(InputEvent.BUTTON4_MASK);
        case 5
            mouse.mousePress(InputEvent.BUTTON5_MASK);
        otherwise
            error("Invalid mouse button");
    end
end

function release(mouse,button)
% Releases specified mouse button.
    import java.awt.Robot;
    % import java.awt.event.*;
    import java.awt.event.KeyEvent;
    import java.awt.event.InputEvent;
    switch button
        case 1
            mouse.mouseRelease(InputEvent.BUTTON1_MASK);
        case 2
            mouse.mouseRelease(InputEvent.BUTTON2_MASK);
        case 3
            mouse.mouseRelease(InputEvent.BUTTON3_MASK);
        case 4
            mouse.mouseRelease(InputEvent.BUTTON4_MASK);
        case 5
            mouse.mouseRelease(InputEvent.BUTTON5_MASK);
        otherwise
            error("Invalid mouse button");
    end
end

function type_keys(mouse,string)
% Types keys of `string` in sequence using Java Robot.
    global timer_mult
    string = char(string);
    uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    for key = string
        % Exception type 1: Character doesn't match KeyEvent
        if isequal(key," ")
            key = "SPACE";      % In case of space in `pid`
%         elseif isequal(key,":")
%             key = "COLON";
        elseif isequal(key,"/")
            key = "SLASH";
        elseif isequal(key,"\")
            key = "BACK_SLASH";
%         elseif isequal(key,"@")
%             key = "AT";
        elseif isequal(key,".")
            key = "PERIOD";
        elseif isequal(key,",")
            key = "COMMA";
        elseif isequal(key,"+")
            key = "PLUS";
        elseif isequal(key,"-")
            key = "MINUS";
%         elseif isequal(key,"_")
%             key = "UNDERSCORE";
        elseif isequal(key,"=")
            key = "EQUALS";
%         elseif isequal(key,"!")
%             key = "EXCLAMATION_MARK";
        end
        % Exception type 2: KeyEvent does not exist
        import java.awt.event.KeyEvent;
        if isequal(key,":")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_SEMICOLON)");
                eval("mouse.keyRelease(KeyEvent.VK_SEMICOLON)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        elseif isequal(key,"@")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_2)");
                eval("mouse.keyRelease(KeyEvent.VK_2)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        elseif isequal(key,"&")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_7)");
                eval("mouse.keyRelease(KeyEvent.VK_7)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        elseif isequal(key,"_")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_MINUS)");
                eval("mouse.keyRelease(KeyEvent.VK_MINUS)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        elseif isequal(key,"!")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_1)");
                eval("mouse.keyRelease(KeyEvent.VK_1)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        elseif isequal(key,"?")
            eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_SLASH)");
                eval("mouse.keyRelease(KeyEvent.VK_SLASH)");
            eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
        else
            % Type character!
            if ismember(key,uppercase)
                % Capitalize the letter, if need be
                eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                    eval("mouse.keyPress(KeyEvent.VK_"+upper(key)+")");
                    eval("mouse.keyRelease(KeyEvent.VK_"+upper(key)+")");
                eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
            else
                % Type character normally
                eval("mouse.keyPress(KeyEvent.VK_SHIFT)");
                eval("mouse.keyRelease(KeyEvent.VK_SHIFT)");
                eval("mouse.keyPress(KeyEvent.VK_"+upper(key)+")");
                eval("mouse.keyRelease(KeyEvent.VK_"+upper(key)+")");
            end
        end
        pause(timer_mult*0.001);
    end
end

% CMD Functions
function cmd_open(varargin)
% Begins running app if it is not already open.
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            after_pause = 0.5;
        case 2
            app_name = varargin{1};
            after_pause = varargin{2};
        otherwise
            error("cmd_open() : Incorrect number of input arguments.");
    end
    path = pwd;
    
    % Is app open?
    isopen = cmd_isopen(app_name);
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'RenderDoc';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'Chromium';    % Also applies to Chromium.
    end
    
    if ~isopen
        command = "start "+path+"\"+app_name;  % Not case-sensitive
        [status,cmdout] = system(command);
        pause(timer_mult*after_pause);
    else
%         % Switch to app
%         cmd_switchto(app_name,0);
        % Do nothing
    end
end

function isopen = cmd_isopen(varargin)
% Checks if app is open.
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            duration = 0.5;
        case 2
            app_name = varargin{1};
            duration = varargin{2};
        otherwise
            error("cmd_wait_until_open() : Incorrect number of input arguments.");
    end
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'notepad')
        app_name = 'notepad';
    elseif isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'qrenderdoc';
        app_name = 'renderdoccmd';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'chrome';    % Also applies to Chromium.
    elseif isequal(lower(app_name),'chromium')
        app_name = 'chrome';    % Also applies to Chrome.
    end
    
    command = append('is_running.bat ',app_name,'.exe');  % YES, case-sensitive! (Must be EXACT match.)
    [status,cmdout] = system(command);
    if contains(cmdout,'Running!!')
        isopen = true;
    else
        isopen = false;
    end
end

function cmd_assert_open(varargin)
% Checks if app is open.
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            duration = 0.5;
        case 2
            app_name = varargin{1};
            duration = varargin{2};
        otherwise
            error("cmd_wait_until_open() : Incorrect number of input arguments.");
    end
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'notepad')
        app_name = 'notepad';
    elseif isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'qrenderdoc';
        app_name = 'renderdoccmd';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'chrome';    % Also applies to Chromium.
    elseif isequal(lower(app_name),'chromium')
        app_name = 'chrome';    % Also applies to Chrome.
    end
    
    command = append('is_running.bat ',app_name,'.exe');  % YES, case-sensitive! (Must be EXACT match.)
    [status,cmdout] = system(command);
    if contains(cmdout,'Running!!')
        isopen = true;
    else
        error(app_name+" is not open.");
    end
end

function cmd_wait_until_open(varargin)
% Waits until app of specified name (without '.exe', case-sensitive) opens,
% or until timeout limit is exceeded (in which case an error is raised).
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            after_pause = 2;
        case 2
            app_name = varargin{1};
            after_pause = varargin{2};
        otherwise
            error("cmd_wait_until_open() : Incorrect number of input arguments.");
    end
    
    % Set timeout duration
    timeout_duration = 30;  % In seconds
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'notepad')
        app_name = 'notepad';
    elseif isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'qrenderdoc';
        app_name = 'renderdoccmd';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'chrome';    % Also applies to Chromium.
    elseif isequal(lower(app_name),'chromium')
        app_name = 'chrome';    % Also applies to Chrome.
    end
    
    tic
    while toc < timeout_duration
        command = append('is_running.bat ',app_name,'.exe');  % YES, case-sensitive! (Must be EXACT match.)
        [status,cmdout] = system(command);
        if contains(cmdout,'Running!!')
            pause(timer_mult*after_pause);
%             disp("RUNNING");
            return
        end
        pause(timer_mult*1);
        disp("...");
    end
    % If timeout duration is exceeded without app opening
    error("cmd_wait_until_open() : '"+app_name+".exe' did not open within timeout limit of "+timeout_duration+" seconds.");
    
end

function cmd_wait_for_response(varargin)
% Waits until app of specified name (without '.exe', case-sensitive)
% responds (at least 6 seconds).
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            after_pause = 0.2;
        case 2
            app_name = varargin{1};
            after_pause = varargin{2};
        otherwise
            error("cmd_wait_until_open() : Incorrect number of input arguments.");
    end
    
    % Set timeout duration
    timeout_duration = 999;  % In seconds
    is_running = cmd_isopen(app_name,0);
    assert(is_running);
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'notepad')
        app_name = 'notepad';
    elseif isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'qrenderdoc';
        app_name = 'renderdoccmd';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'chrome';    % Also applies to Chromium.
    elseif isequal(lower(app_name),'chromium')
        app_name = 'chrome';    % Also applies to Chrome.
    end
    
    pause(timer_mult*6); % Windows takes ~5 seconds for an application to "not respond".
    tic
    while toc < timeout_duration
        command = append('is_responding.bat ',app_name,'.exe');  % YES, case-sensitive! (Must be EXACT match.)
        [status,cmdout] = system(command);
        if ~contains(cmdout,"Application not responding!!")
            pause(timer_mult*after_pause);
            disp("Ready!");
            return
        end
        pause(timer_mult*1);
        disp("... waiting for response from '"+app_name+"' ...");
    end
    % If timeout duration is exceeded without app opening
    error("cmd_wait_for_response() : '"+app_name+".exe' did not open within timeout limit of "+timeout_duration+" seconds.");
    
end

function cmd_switchto(varargin)
% Waits until app of specified name (without '.exe', case-sensitive) opens,
% and then switches to it.
% If unsuccessful, retries repeatedly until timeout limit is exceeded.
    global timer_mult
    switch nargin
        case 1
            app_name = varargin{1};
            after_pause = 0.2;
        case 2
            app_name = varargin{1};
            after_pause = varargin{2};
        otherwise
            error("cmd_switchto() : Incorrect number of input arguments.");
    end
    
    % Substitute commonly used app names
    if isequal(lower(app_name),'notepad')
        app_name = 'notepad';
    elseif isequal(lower(app_name),'blender')
        app_name = 'blender';
    elseif isequal(lower(app_name),'renderdoc')
        app_name = 'renderdoc';
    elseif isequal(lower(app_name),'chrome')
        app_name = 'chrome';    % Also applies to Chromium.
    end
    
    
%     cmd_wait_until_open(app_name,0.2);
    cmd_assert_open(app_name,0.2);
    command = append('switch.vbs "',app_name,'"');  % NOT case-sensitive!
    [status,cmdout] = system(command);
    pause(timer_mult*after_pause);
    
    % Tap: Win+Up (Maximize)
    % Initialize Java Robot
    import java.awt.Robot;
    mouse = Robot;
    % import java.awt.event.*;
    import java.awt.event.KeyEvent;
    mouse.keyPress(KeyEvent.VK_WINDOWS)
        pause(timer_mult*0.1)
        mouse.keyPress(KeyEvent.VK_UP)
        mouse.keyRelease(KeyEvent.VK_UP)
    mouse.keyRelease(KeyEvent.VK_WINDOWS)
    pause(timer_mult*0.1);
end