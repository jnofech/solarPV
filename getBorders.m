function [rbutt_x,rbutt_y,varargout] = getBorders(global_left, global_top, global_width, global_height,message)
% A modified version of `getButton.m` from CodingLikeMad Youtube Channel, 2019

    % Get the button coordinates.
    findWindow = false;

    [imgData] = takeScreenshot(global_left, global_top, global_width, global_height);
    imshow(imgData)
    title(message);

    while findWindow == false
        [rbutt_x,rbutt_y] = ginput(2);
        rbuttCoorX = (rbutt_x(2)+rbutt_x(1))/2;
        rbuttCoorY = (rbutt_y(2)+rbutt_y(1))/2;
        [imgData] = imgData(rbutt_y(1):rbutt_y(2), rbutt_x(1):rbutt_x(2),:);
        imshow(imgData)
        result = questdlg('Was this correct?', 'Instructions', 'No');
        switch result
            case 'Yes'
                findWindow = true;
            case 'No'
                [imgData] = takeScreenshot(global_left, global_top, global_width, global_height);
                imshow(imgData)
            case 'Cancel'
                error('Ending execution.');
        end            
    end
    if nargout==3
        varargout{1} = imgData;
    end
end