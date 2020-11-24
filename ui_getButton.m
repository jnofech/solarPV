function [rbuttCoorX,rbuttCoorY] = ui_getButton(varargin)
% Returns user-inputted coordinates of some point on a screenshot of the
% entire screen.
    switch nargin
        case 0
            message = " ";
            snippet = [];
            ratio = 6;  % Ratio of "screenshot" to "snippet" in figure.
        case 1
            message = varargin{1};
            snippet = [];
            ratio = 6;  % Ratio of "screenshot" to "snippet" in figure.
        case 2
            message = varargin{1};
            snippet = varargin{2};
            ratio = 6;  % Ratio of "screenshot" to "snippet" in figure.
        case 3
            message = varargin{1};
            snippet = varargin{2};
            ratio   = varargin{3};  % Ratio of "screenshot" to "snippet" in figure.
        otherwise
            error("ui_getButton() : Too many arguments.")
    end


    % Get screen dimensions
    screenSize = get(0, 'screensize');
    left = screenSize(1);
    top  = screenSize(2);
    width  = screenSize(3);
    height = screenSize(4);
    
    button_found = false;
    while button_found==false
        % Take Screenshot
        img = takeScreenshot(left,top,width,height);

        % Display screenshot
        figure(1);
        if ~isempty(snippet)
            subplot(ratio,1,1);
            imshow(snippet);
            subplot(ratio,1,2:ratio);
        end
        imshow(img);
        set(gcf,'position',[left,top,width,height]);

        % Input text box location
        title(message);
        [rbuttCoorX,rbuttCoorY] = ginput(1);
        rbuttCoorX = round(rbuttCoorX);
        rbuttCoorY = round(rbuttCoorY);
        close all;

        % Display inputted result
        img_snippet = img(max(rbuttCoorY-25,top):min(rbuttCoorY+25,height),max(rbuttCoorX-100,left):min(rbuttCoorX+100,width),:);
        imshow(img_snippet); axis on;
        hold on;
        rbuttCoorY_plot = size(img_snippet,1) * (rbuttCoorY - max(rbuttCoorY-25,top)) / (min(rbuttCoorY+25,height) - max(rbuttCoorY-25,top));
        rbuttCoorX_plot = size(img_snippet,2) * (rbuttCoorX - max(rbuttCoorX-100,left)) / (min(rbuttCoorX+100,width) - max(rbuttCoorX-100,left));
        plot(rbuttCoorX_plot, rbuttCoorY_plot,'rx','MarkerSize',10);
        hold off;

        % Confirm usability!
        result = questdlg('Is this correct?', 'Confirmation','Yes', 'No', 'Cancel', 'Yes');
        switch result
            case 'Yes'
                button_found = true;
            case 'No'
                % Do nothing; i.e. retry
            case 'Cancel'
                error('Ending execution.');
        end
        close all;
    end
end