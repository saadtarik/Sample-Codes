%ENCE689G PS#6 (Problem 1)
%
%Author: Saad Tarik
%Created: Oct 24, 2014

clear;clc;tic

start_dir = pwd;
load state_outlines_Matlab.mat; 
% load coast;

%% Specify user-defined/graphical parameters
image_type = '-depsc';
custom_colors = [1 1 1; jet(100)];
afs = 10;
state_outlines = shaperead('state_outlines','UseGeoCoords',true);
land_outlines = shaperead('landareas','UseGeoCoords',true);

%% Part (a)
% Examine the composite reflectivity file
filename = 'KLWX_NCR_20140812_000400.nc';
% filename = 'KLWX_DPR_20140812_000400.nc';
% filename = 'KLWX_NHL_20140812_000400.nc';
% filename = 'KLWX_NML_20140812_000400.nc';
% filename = 'KLWX_NLL_20140812_000400.nc';
ncdisp(filename)

%% Part (b)
cref = ncread('KLWX_NCR_20140812_000400.nc','cref');
lon = ncread('KLWX_NCR_20140812_000400.nc','lon');
lat = ncread('KLWX_NCR_20140812_000400.nc','lat');
[lon2D,lat2D] = meshgrid(lon,lat);
precip = ncread('KLWX_DPR_20140812_000400.nc','dpr');

% Start graphical output
figure(1);clf
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,cref','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)

subplot(1,2,2)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,precip','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Instantaneous Precipitation Rate [in/hr]')
colorbar
set(gcf,'Colormap',custom_colors)

%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r300','PS06_Prob_1b');close

%% Part (c)
hcref = ncread('KLWX_NHL_20140812_000400.nc','hcref');
mcref = ncread('KLWX_NML_20140812_000400.nc','mcref');
lcref = ncread('KLWX_NLL_20140812_000400.nc','lcref');

figure(2);clf
set(gcf,'Position',[100 100 800 600])
subplot(2,2,1)
% ax = usamap({'MD','VA','DC'});
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,hcref','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('High-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)
subplot(2,2,2)
% ax = usamap({'MD','VA','DC'});
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,mcref','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Mid-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)
subplot(2,2,3)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,lcref','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Low-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r300','PS06_Prob_1c');close

%% Part (d)
cref2 = ncread('KLWX_NCR_20140812_203400.nc','cref');
precip2 = ncread('KLWX_DPR_20140812_203400.nc','dpr');
hcref2 = ncread('KLWX_NHL_20140812_203400.nc','hcref');
mcref2 = ncread('KLWX_NML_20140812_203400.nc','mcref');
lcref2 = ncread('KLWX_NLL_20140812_203400.nc','lcref');

figure(3);clf
set(gcf,'Position',[100 100 1000 800])
subplot(3,2,1)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,cref2','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)

subplot(3,2,2)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,precip2','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Instantaneous Precipitation Rate [in/hr]')
colorbar
set(gcf,'Colormap',custom_colors)

subplot(3,2,3)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,hcref2','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('High-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)

subplot(3,2,4)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,mcref2','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Mid-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)

subplot(3,2,5)
ax = worldmap([36 42],[-80 -74]);
geoshow(lat2D,lon2D,lcref2','DisplayType','surface')
% Display state outlines
geoshow([state_outlines.Lat],[state_outlines.Lon],'Color','black',...
    'LineWidth',0.2)
title('Low-Level Composite Reflectivity [dBZ]')
caxis([0 40]);colorbar
set(gcf,'Colormap',custom_colors)
%save figure for later use
set(gcf,'PaperPositionMode','auto')
print(image_type,'-r300','PS06_Prob_1d');close

%% Echo computational run-time
fct_echo_runtime('seconds')
