% ENCE689G PS#8 (Prob 2)
%
% Author: Saad Tarik
% Created: Nov 12, 2014

clear; clc; tic
start_dir = pwd;
load Northern_Hemisphere_land_sea_ice_mask.mat;
load MERRA_rainfall_20090801_to_20090803.mat;
load NL_geolocation_arrays.mat;
coast = load('coast');
% Preprocess land_sea_ice_mask
ind0 = find(land_sea_ice_mask == 0);
land_sea_ice_mask(ind0) = 1;
indnan = find(land_sea_ice_mask ~= 1);
land_sea_ice_mask(indnan) = NaN;

%% Specify user-defined and graphical paramteres
image_type = '-depsc';
image_res = '-r200';
font_size = 11;
custom_colors = [1 1 1;jet(500)];

%% Part (a)
filename = 'ID2r3-AMSRE-NL2009213A.v03.10H';
Tb_10H = fct_load_AMSR_E_Tb(filename);
Tb_10H_land = Tb_10H .* land_sea_ice_mask;
filename = 'ID2r3-AMSRE-NL2009213A.v03.10V';
Tb_10V = fct_load_AMSR_E_Tb(filename);
Tb_10V_land = Tb_10V .* land_sea_ice_mask;
PR_213 = (Tb_10V_land - Tb_10H_land) ./ (Tb_10V_land + Tb_10H_land);

% Start Visualization
figure(1);clf
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,Tb_10H_land,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([220 300])
title('H-pol','FontSize',font_size)
title(h,'Tb [K]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)

subplot(1,2,2)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,Tb_10V_land,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([220 300])
title('V-pol','FontSize',font_size)
title(h,'Tb [K]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_2a'),close

%% Part (b)
filename = 'ID2r3-AMSRE-NL2009215A.v03.10H';
Tb_10H = fct_load_AMSR_E_Tb(filename);
Tb_10H_land = Tb_10H .* land_sea_ice_mask;
filename = 'ID2r3-AMSRE-NL2009215A.v03.10V';
Tb_10V = fct_load_AMSR_E_Tb(filename);
Tb_10V_land = Tb_10V .* land_sea_ice_mask;
PR_215 = (Tb_10V_land - Tb_10H_land) ./ (Tb_10V_land + Tb_10H_land);

% Start Visualization
figure(2);clf
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,Tb_10H_land,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([220 300])
title('H-pol','FontSize',font_size)
title(h,'Tb [K]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)

subplot(1,2,2)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,Tb_10V_land,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([220 300])
title('V-pol','FontSize',font_size)
title(h,'Tb [K]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_2b'),close

%% Part (c)
figure(3);clf
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,PR_213,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([-0.1 0.1])
title('Day 213','FontSize',font_size)
title(h,'PR [-]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)

subplot(1,2,2)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,PR_215,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([-0.1 0.1])
title('Day 215','FontSize',font_size)
title(h,'PR [-]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_2c'),close

%% Part (d)
dPR = PR_215 - PR_213;
rainfall_cm = rainfall_cm .* land_sea_ice_mask;

% Start Visualization
figure(4);clf
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,dPR,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([-0.025 0.025])
title('\DeltaPR','FontSize',font_size)
title(h,'\DeltaPR [-]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)

subplot(1,2,2)
axesm('stereo','MapLatLimit',[0 90]);
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-30)
geoshow(EASE_Grid.lat2D,EASE_Grid.lon2D,rainfall_cm,'DisplayType','surface')
hold on
plotm(coast.lat,coast.long,'-k')
h = colorbar;caxis([0 2])
title('Cumulative Rainfall','FontSize',font_size)
title(h,'[cm]','FontSize',font_size)
set(gcf,'ColorMap',custom_colors)
%save image to file
set(gcf,'PaperPositionMode','auto')
print(image_type,image_res,'PS08_Prob_2d'),close

%% Echo computational run-time
fct_echo_runtime('seconds')
