% ENCE689G PS#8 (Prob 3)
%
% Author: Saad Tarik
% Created: Nov 12, 2014

clear; clc; tic
start_dir = pwd;
%load Northern_Hemisphere_land_sea_ice_mask.mat;
%load MERRA_rainfall_20090801_to_20090803.mat;
%load NL_geolocation_arrays.mat;
coast = load('coast');
% Preprocess land_sea_ice_mask
% ind0 = find(land_sea_ice_mask == 0);
% land_sea_ice_mask(ind0) = 1;
% indnan = find(land_sea_ice_mask ~= 1);
% land_sea_ice_mask(indnan) = NaN;

%% Specify user-defined and graphical paramteres
image_type = '-depsc';
image_res = '-r200';
font_size = 11;
custom_colors = [1 1 1;jet(500)];

%% Part (a)
filename = 'GRCTellus.JPL.200204_201406.LND.RL05.DSTvSCS1401.nc';
lwe_thickness = ncread(filename,'lwe_thickness');
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
time_vec = ncread(filename,'time');
lon = circshift(lon,[0 180]); % Shift to match w/ panoply plot
lwe_t1 = lwe_thickness(:,:,1);
lwe_tend = lwe_thickness(:,:,136);
[lon2D,lat2D] = meshgrid(lon,lat);

% Start Visualization
figure(1);clf
set(gcf,'Position',[100 100 1200 600])
subplot(1,2,1)
ax = worldmap([-90 90],[-180 180]);
geoshow(lat2D,lon2D,lwe_t1','DisplayType','surface')
title('2 April 2002','FontSize',font_size)
h_cb = colorbar('southoutside');title(h_cb,'Liquid Water Equivalent [cm]')
caxis([-100 100]);set(gcf,'ColorMap',custom_colors); gridm on
hold on
plotm(coast.lat,coast.long,'-k')

subplot(1,2,2)
ax = worldmap([-90 90],[-180 180]);
geoshow(lat2D,lon2D,lwe_tend','DisplayType','surface')
title('2 June 2014','FontSize',font_size)
h_cb = colorbar('southoutside');title(h_cb,'Liquid Water Equivalent [cm]')
caxis([-100 100]);set(gcf,'ColorMap',custom_colors); gridm on
hold on
plotm(coast.lat,coast.long,'-k')

%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_3b'),close

%% Part (c)
st_dev = std(lwe_thickness,0,3);

% Start Visualization
figure(2);clf
set(gcf,'Position',[100 100 600 400])
ax = worldmap([-90 90],[-180 180]);
geoshow(lat2D,lon2D,st_dev','DisplayType','surface')
h_cb = colorbar('southoutside');title(h_cb,'Standard Deviation [cm]')
caxis([0 30]);set(gcf,'ColorMap',custom_colors); gridm on
hold on
plotm(coast.lat,coast.long,'-k')

%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_3c'),close

%% Part (d)
lwe_Amazon = lwe_thickness(300,85,:);
date_vector = time_vec + datenum(2002,1,1) - 1;

% Start Visualization
figure(3);clf
plot(date_vector,squeeze(lwe_Amazon),'-k')
datetick('x',1)
xlabel('Date');ylabel('Liquid Water Equivalent [cm]')

%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_3d'),close

%% Part (e)
long_ave_lwe = squeeze(nanmean(lwe_thickness,1));

% Start visualization
figure(4);clf
imagesc(date_vector,lat,long_ave_lwe)
set(gca,'YDir','normal');colorbar
datetick('x',1)
xlim([date_vector(1) date_vector(end)])
xlabel('Date');ylabel('Latitude [deg]')

%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_3e'),close

%% Echo computational run-time
fct_echo_runtime('seconds')
