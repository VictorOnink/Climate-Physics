% Code for plotting just the geostrophic velocities, not the eddies on top
% of them

close all
clear all
% see netCDF file description
%ncdisp('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995')

% load track data
load('AllTracks.mat')
time=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','time');
vgos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','vgos');
ugos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','ugos');
aVgos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','vgosa');
aUgos=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','ugosa');
lat=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','latitude');
lon=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/GeostrophicVelocityMozam1993_1995','longitude');
[Lat, Lon]=meshgrid(lat, lon);
%Where the images will be saved
workingDir = '/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies';
mkdir(workingDir)
%Create unit vectors
magnitude=sqrt(aVgos.^2+aUgos.^2);
aVnorm=(aVgos./magnitude);
aUnorm=(aUgos./magnitude);
Vorticity=zeros(length(lat),length(lon),length(time));
for i=1:length(time)
   %Compute the vorticity at each point
   [aspectU, slopeU,gradNU,gradEU]=gradientm(Lat,Lon,ugos(:,:,i));
   [aspectV, slopeV,gradNV,gradEV]=gradientm(Lat,Lon,vgos(:,:,i));
   Vorticity(:,:,i)=gradEV-gradNU;
   %[Vorticity(:,:,i),velocity]=curl(Lon,Lat,ugos(:,:,i),vgos(:,:,i));
end
plusVorU=aUnorm;
plusVorU(Vorticity<0)=NaN;
plusVorV=aVnorm;
plusVorV(Vorticity<0)=NaN;
negVorU=aUnorm;
negVorU(Vorticity>0)=NaN;
negVorV=aVnorm;
negVorV(Vorticity>0)=NaN;
%%
%Define the total domain of the Mozambique Channel
yv_tot=[-10,-10,-12,-25.6,-26.5,-10];   %Degrees latitude south
xv_tot=[32,40, 49.2, 45, 32,32];    %Degrees longitude East


%%
close all
%set the background map of the animation
j=1;
%folder into which we will save the images
folder = '/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies';

for i=1%length(time)
    tic
    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    map=axesm('miller','maplatlim',[-30, -10],'maplonlim',[30, 50],'grid','on',...
        'meridianlabel','on','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','on','plabellocation',[-10 -15 -20 -25 -30 -35],...
        'fontsize',12);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green'); 
    title({'Geostrophic Velocity Anomalies Mozambique Channel';'01-01-1993 to 01-01-1995'},'fontsize',18)
    [year, month, day]=julian2greg(time(i)+2433282.5);
    y=get(gca,'ylim');
    x=get(gca,'xlim');
    string=sprintf('Date: %d-%d-%d',day,month,year);
    dim = [.325 .62 .3 .3];
    annotation('textbox',dim,'String',string,'FitBoxToText','on','FontSize',16,...
        'BackgroundColor','w');
    
    quivermc(Lat,Lon,plusVorU(:,:,i),plusVorV(:,:,i),'color','r');
    quivermc(Lat,Lon,negVorU(:,:,i),negVorV(:,:,i),'color','b');
	% Create a PNG filename.
	pngFileName = sprintf('MozambiqueChannelAnom19931995_%d.png', i);
	fullFileName = fullfile(folder, pngFileName);
	saveas(gcf, fullFileName);
    %close all
    toc
end

%% Create a video of the geostrophic velocities
a=dir(['/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies' '/*.png']);
frames=cell(size(a,1),1);
%Create working directory
workingDir = '/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies/Movie';
mkdir(workingDir)
%create AVI file with all the figures
writerObj = VideoWriter('/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies/Movie/GeostrophicAnomEddy19931995.avi');
writerObj.FrameRate = 15;
open(writerObj);
for i = 1 : size(a,1)
  GeoThatDay = imread(sprintf('/Users/victoronink/Desktop/DataFiles/MAIO/JustVelocityPlots/Anomalies/MozambiqueChannelAnom19931995_%d.png', i));
  writeVideo(writerObj, GeoThatDay);
end
close(writerObj);