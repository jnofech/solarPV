%% TO-DO TOMORROW:
%   * Continue implementing ROI functionality into an app interface. Learn
%     how to input parameters and output results!

%% ROI test

hmap = b.Z;

imagesc(hmap); axis image; colorbar


% Create plot
xmax_plot = size(hmap,2);    ymax_plot = size(hmap,1);
width=1366;height=768;
set(gcf,'position',[200,200,width,height]);
imagesc(hmap, 'XData', (0:(xmax_plot-1))*b.pix2m, 'YData', (0:(ymax_plot-1))*b.pix2m);
colorbar;
axis image;
xlabel('X (m)'); ylabel('Y (m)'); colormap(brewermap([],'RdYlBu'));

% ROI input
roi = drawellipse('StripeColor','y');

%% Peppers test

rgb = imread('peppers.png');
imshow(rgb);
I = rgb2gray(rgb);
hold on
h = imshow(I); % Save the handle; we'll need it later
hold off

[rows,cols] = size(I);
block_size = 200;
P = ceil(rows / block_size);
Q = ceil(cols / block_size);

alpha_data = checkerboard(block_size, P, Q) > 0;
alpha_data = alpha_data(1:rows, 1:cols);
set(h, 'AlphaData', alpha_data);


%% Test - Concatenated Colourmaps
hmap = b.Z;
hmap = hmap - min(hmap,[],'all');   % Lowest value must be zero.

tic
% Colormap selection
cmap_length = 64;
cmap_good = brewermap([cmap_length],'RdYlBu');  % Selected
cmap_bad = gray(cmap_length);                   % Deselected
cmap_combined = [cmap_good; cmap_bad];          % Displays both together

% Plot parameters
xmax_plot = size(hmap,2);    ymax_plot = size(hmap,1);
width=1366;height=768;

% Draw initial map
set(gcf,'position',[200,200,width,height]);
imshow(hmap/max(hmap,[],'all'),'colormap',cmap_good);
rect = drawrectangle();
% rect = drawrectangle('Position',[1,1,1,1]);
deselected = ~createMask(rect);      % Binary of deselected regions. To be determined via app
% deselected = zeros(size(hmap));
% deselected(400:1000,400:1200) = 1;

% % Create "workaround" heightmap
% hmap_disp = hmap;
% hmap_disp = hmap_disp + deselected*max(hmap,[],'all');
% hmap_disp(~isfinite(hmap_disp)) = max(hmap,[],'all');
% 
% % Create plot
% % imagesc(hmap_disp, 'XData', (0:(xmax_plot-1))*b.pix2m, 'YData', (0:(ymax_plot-1))*b.pix2m);
% % imshow(hmap_disp/(2*max(hmap,[],'all'));
% imshow(hmap_disp/(2*max(hmap,[],'all')),'colormap',cmap_combined);
% xlabel('X (m)'); ylabel('Y (m)'); 
% axis image;

toc

%% Test function using actual app!
% Parameters
name = b.name;
hmap = b.Z;
map = b.map_r;

hmap = hmap - min(hmap,[],'all');   % Lowest value must be zero.
outputname = "output/"+name+"_"+pix_length+"pix_mask.mat";
if ~isfile(outputname)
    % Create output mask!
    roi_app_test(name,hmap,map,outputname);
end
load(outputname);


%% Test the app-test output!
load(outputname);
% Create "workaround" heightmap
hmap_disp = hmap;
hmap_disp = hmap_disp + ~mask*max(hmap,[],'all');
hmap_disp(~isfinite(hmap_disp)) = max(hmap,[],'all');

% Create plot
% imagesc(hmap_disp, 'XData', (0:(xmax_plot-1))*b.pix2m, 'YData', (0:(ymax_plot-1))*b.pix2m);
% imshow(hmap_disp/(2*max(hmap,[],'all'));
imshow(hmap_disp/(2*max(hmap,[],'all')),'colormap',cmap_combined); colorbar;
xlabel('X (m)'); ylabel('Y (m)'); 
axis image;