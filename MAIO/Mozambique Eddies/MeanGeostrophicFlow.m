% Create plot of the average flow speed over the period of 1991-1993
clear all
% define below the coordinates of your framework
    lon_min = 30; lon_max = 50;
    lat_min  = -30;  lat_max = -10;
    
% extract bathymetry for the region of study    
[elevation, vlat, vlon] = mygrid_sand_1m([lat_min lat_max lon_min lon_max],1);

%Data
time=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','time');
vgos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','vgos');
ugos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','ugos');
lat=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','latitude');
lon=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','longitude');
[Lat, Lon]=meshgrid(lat(1:2:end), lon(1:2:end));
%Means over 2 year period
meanVgos=mean(vgos,3);
meanUgos=mean(ugos,3);
magnitudes=sqrt(meanVgos.^2+meanUgos.^2);
%%
close all
figure(1)
title('Mean Geostrophic Flow 01-01-1993 to 01-01-1995','FontSize',16)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% initiate a map projection
    m_proj('mercator','lat',[lat_min lat_max],'lon',[lon_min lon_max]); hold on    
% plot bathymetry
    [c1,h1] = m_contourf(vlon,vlat,elevation,[-1000 -2000 -3000 -4000 -5000 -6000]); hold on
    set(h1,'color','k');
    colormap('gray')
    c=colorbar;
    c.Label.String = 'Depth contours (m)';
    c.FontSize= 12;
% call the coastline previously saved
    m_usercoast('/Users/victoronink/Desktop/DataFiles/MAIO/moz_coast.mat','patch',[.7 .7 .7]);
    m_grid('linestyle','none','tickedir','out','fontsize',14);
% Plot vector grid of mean flow velocity
    keystr='m/s';
    Patch=m_vec(1,Lon,Lat,meanUgos(1:2:end,1:2:end),meanVgos(1:2:end,1:2:end),...
        'r','shaftwidth',0.2,'edgeclip','on');
    
folder = '/Users/victoronink/Desktop/DataFiles/MAIO';
fullFileName = fullfile(folder, 'MeanGeoFlow19931995.png');
saveas(gcf, fullFileName);
