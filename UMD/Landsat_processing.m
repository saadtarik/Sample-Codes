%ENCE689G PS#6 (Problem 3)
%
%Author: Saad Tarik
%Created: Oct 22, 2014

clear;clc;tic
start_dir = pwd;

%% Specify user-defined inputs
dir = 'PS06_Problem02/LE70150332014186EDC00_2/';

%% Specify graphical parameters
image_type = '-depsc';
custom_colors = [1 1 1; jet(100)];
afs = 10;

eval(['cd ',dir])

%% Part(a)

filename = 'LE70150332014186EDC00_B1.TIF';
[Band_1,R] = geotiffread(filename);
Band_1 = double(Band_1);
ind_NaN = find(Band_1 == 0);
Band_1(ind_NaN) = NaN;
[nrows, ncols] = size(Band_1);
easting = linspace(R.XWorldLimits(1),R.XWorldLimits(2),nrows);
northing = linspace(R.YWorldLimits(1),R.YWorldLimits(2),ncols);

filename = 'LE70150332014186EDC00_B3.TIF';
[Band_3,R] = geotiffread(filename);
Band_3 = double(Band_3);
ind_NaN = find(Band_3 == 0);
Band_3(ind_NaN) = NaN;

filename = 'LE70150332014186EDC00_B4.TIF';
[Band_4,R] = geotiffread(filename);
Band_4 = double(Band_4);
ind_NaN = find(Band_4 == 0);
Band_4(ind_NaN) = NaN;

% Start graphical output
figure(1);clf
subplot(2,2,1)
imagesc(easting,northing,Band_1)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Landsat (Band 1)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,2)
imagesc(easting,northing,Band_3)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Landsat (Band 3)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,3)
imagesc(easting,northing,Band_4)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Landsat (Band 4)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)

%save figure for later use
%set(gcf,'PaperPositionMode','auto')
%print(image_type,'-r100','PS06_Prob_3a');close

%% Prob 3(b)
LL_B1 = -7.38071 + Band_1 * 1.181;
LL_B3 = -5.94252 + Band_3 * 0.943;
LL_B4 = -6.06929 + Band_4 * 0.969;

% Start graphical output
figure(2);clf
subplot(2,2,1)
imagesc(easting,northing,LL_B1)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Radiance (Band 1)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,2)
imagesc(easting,northing,LL_B3)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Radiance (Band 3)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,3)
imagesc(easting,northing,LL_B4)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Radiance (Band 4)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r100','PS06_Prob_3b');close

%% Prob 3(c)
theta = 0.2898; %[rad]
d_ES = 1.01667; %[astronomical units]
ESUN_B1 = 1997;
ESUN_B3 = 1533;
ESUN_B4 = 1039;

r_B1 = (pi * LL_B1 * d_ES^2) / (ESUN_B1 * cos(theta));
r_B3 = (pi * LL_B3 * d_ES^2) / (ESUN_B3 * cos(theta));
r_B4 = (pi * LL_B4 * d_ES^2) / (ESUN_B4 * cos(theta));

% Start graphical output
figure(3);clf
subplot(2,2,1)
imagesc(easting,northing,r_B1)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Reflectance (Band 1)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,2)
imagesc(easting,northing,r_B3)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Reflectance (Band 3)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
subplot(2,2,3)
imagesc(easting,northing,r_B4)
xlabel('Easting [m]');ylabel('Northing [m]')
title('Reflectance (Band 4)')
colorbar;axis square
set(gcf,'Colormap',custom_colors)
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r100','PS06_Prob_3c');close

%% Prob 3(d)
NDVI = (r_B4-r_B3)./(r_B4+r_B3);
custom_colors_NDVI = [1 1 1; winter(200)];
%Start graphical output
figure(4);clf
imagesc(easting,northing,NDVI)
xlabel('Easting [m]');ylabel('Northing [m]')
title('NDVI [-]')
colorbar;axis square;caxis([-0.9 1])
set(gcf,'Colormap',custom_colors_NDVI)
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r100','PS06_Prob_3d');close

%% Prob 3(e)
%Start graphical output
figure(5);clf
hist(NDVI(:),200)
xlim([-0.9 1])
xlabel('NDVI [-]');ylabel('Number of Occurences')
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r100','PS06_Prob_3e');close
