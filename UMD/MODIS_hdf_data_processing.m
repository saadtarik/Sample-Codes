% ENCE689G PS#7 (Prob 2)
%
% Author: Saad Tarik
% Created: Nov 4, 2014

clear;clc;tic

start_dir = pwd;

%% Specify user-defined and graphical parameters
image_type = '-depsc';
image_res = '-r100';
font_size = 9;
custom_colors = [1 1 1;jet(300)];
custom_colors_bin = [1 1 1;jet(1)];

%% Load file(s) and import variable(s)
load coast;
filename = 'MOD09GA.A2013067.h11v04.005.2013069055208.hdf';
% Part b
dataset = 'sur_refl_b01_1'; ref1 = double(hdfread(filename,dataset));
dataset = 'sur_refl_b02_1'; ref2 = double(hdfread(filename,dataset));
dataset = 'sur_refl_b03_1'; ref3 = double(hdfread(filename,dataset));
dataset = 'sur_refl_b04_1'; ref4 = double(hdfread(filename,dataset));
dataset = 'sur_refl_b06_1'; ref6 = double(hdfread(filename,dataset));

ref1 = (ref1 + 100) ./ 10000;
ref2 = (ref2 + 100) ./ 10000;
ref3 = (ref3 + 100) ./ 10000;
ref4 = (ref4 + 100) ./ 10000;
ref6 = (ref6 + 100) ./ 10000;

true_color(:,:,1) = ref1;true_color(:,:,2) = ref4;true_color(:,:,3) = ref3;

true_color_corrected = true_color;
ind_too_big = find(true_color>1);
true_color_corrected(ind_too_big) = NaN;
ind_too_small = find(true_color<0);
true_color_corrected(ind_too_small) = NaN;

%Start visualization
figure(1);clf
set(gcf,'Position',[100 100 400 400])
imagesc(true_color_corrected)
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS07_Prob_2b'),close

% Part c
NDSI = (ref4 - ref6) ./ (ref4 + ref6);

NDSI_bin = NDSI;
ind_snow = find(NDSI_bin > 0.4);
NDSI_bin(ind_snow) = 1;
ind_no_snow = find(NDSI_bin < 0.4);
NDSI_bin(ind_no_snow) = 0;

%Start visualization
figure(2);clf
set(gcf,'Position',[100 100 1200 400])
subplot(1,3,1)
imagesc(true_color_corrected); axis square
title('True Color')
subplot(1,3,2)
imagesc(NDSI);axis square
title('NDSI')
colorbar;set(gcf,'ColorMap',custom_colors)
subplot(1,3,3)
imagesc(NDSI_bin); axis square
title('NDSI (binary)')
h = colorbar('YTick',0:1:1);set(h,'ColorMap',custom_colors_bin)
set(h,'YTickLabel',{'No Snow','Snow'})
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS07_Prob_2c'),close

% Part d
NDVI = (ref2 - ref1) ./ (ref2 + ref1);

NDVI_bin = NDVI;
ind_veg = find(NDVI_bin > 0.5);
NDVI_bin(ind_veg) = 1;
ind_no_veg = find(NDVI_bin < 0.5);
NDVI_bin(ind_no_veg) = 0;

%Start visualization
figure(3);clf
set(gcf,'Position',[100 100 1200 400])
subplot(1,3,1)
imagesc(true_color_corrected); axis square
title('True Color')
subplot(1,3,2)
imagesc(NDVI);axis square
title('NDVI')
colorbar;set(gcf,'ColorMap',custom_colors)
subplot(1,3,3)
imagesc(NDVI_bin);axis square
title('NDVI (binary)')
h = colorbar('YTick',0:1:1);
set(h,'YTickLabel',{'No Vegetation','Vegetation'})
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS07_Prob_2d'),close
%% Return to original directory and echo computational run-time
% eval(['cd ',start_dir])
fct_echo_runtime('seconds')
