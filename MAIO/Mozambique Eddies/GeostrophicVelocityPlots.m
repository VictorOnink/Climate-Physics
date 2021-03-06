% Code to read in the geostrophic velocities for 01-01-1993 to 01-01-1995
% and to see how they match the tracks of our eddies

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
workingDir = '/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots/Anomalies';
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
end
plusVorU=aUgos;
plusVorU(Vorticity<0)=NaN;
plusVorV=aVgos;
plusVorV(Vorticity<0)=NaN;
negVorU=aUgos;
negVorU(Vorticity>0)=NaN;
negVorV=aVgos;
negVorV(Vorticity>0)=NaN;
%%
%Define the total domain of the Mozambique Channel
yv_tot=[-10,-10,-12,-25.6,-26.5,-10];   %Degrees latitude south
xv_tot=[32,40, 49.2, 45, 32,32];    %Degrees longitude East

%find the tracks that are in the domain in this point
start=15706;
stop=16436;
match={};
m=1;

for i=1:length(Tracks)
    eddy=Tracks{i};
    if eddy(1,8)>=start && eddy(1,8)<=stop
        match{m}=eddy(:,:);
        m=m+1;
    end      
end
% Break separated tracks by lifetime, so only those that last longer than
% 60 days
agesplit=60; %days
Match30plus=[];
count=1;
for i=1:length(match)
    track=match{i};
    if length(track)>=agesplit
        Match30plus{count}=track;
        count=count+1;
    end
end
starttime=[];
for i=1:length(Match30plus)
    track=Match30plus{i};
    starttime(end+1,:)=[track(1,8),length(track(:,8))];
end
%%
close all

j=1;
%folder into which we will save the images
folder = '/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots';
for i=1:length(time)
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
    for j=1:length(Match30plus)
        track=Match30plus{j};
        if time(i)>=starttime(j,1) && time(i)<=(starttime(j,1)+starttime(j,2)-1)
            Eddies = strcat('eddy',num2str(j));
            k=find(time(i)==track(:,8));
            if track(k,2)==1 %anticyclonic red
                EddyCircle.(Eddies)=circlem(track(k,3),track(k,4),track(k,7),'facecolor','r',...
                    'facealpha',.15,'edgecolor','none');
            elseif track(k,2)==-1 %cyclonic blue
                EddyCircle.(Eddies)=circlem(track(k,3),track(k,4),track(k,7),'facecolor','b',...
                    'facealpha',.15,'edgecolor','none');
            end
        end   
    end
    quivermc(Lat,Lon,plusVorU(:,:,i),plusVorV(:,:,i),'color','r','reference',0.8);
    quivermc(Lat,Lon,negVorU(:,:,i),negVorV(:,:,i),'color','b','reference',0.8);
    
	% Create a PNG filename.
	pngFileName = sprintf('MozambiqueEddy19931995NotNorm_%d.png', i);
	fullFileName = fullfile(folder, pngFileName);
	saveas(gcf, fullFileName);
    close all
    toc
end

%% Create a video of the geostrophic velocities with eddies
a=dir(['/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots/Anomalies' '/*.png']);
frames=cell(size(a,1),1);
%Create working directory
workingDir = '/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots/Anomalies/Movie';
mkdir(workingDir)

%create AVI file with all the figures
writerObj = VideoWriter('/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots/Anomalies/Movie/GeostrophicAnomEddiesNotNorm19931995.avi');
writerObj.FrameRate = 5;
open(writerObj);
for i = 1 : size(a,1)
  GeoThatDay = imread(sprintf('/Users/victoronink/Desktop/DataFiles/MAIO/VelocityPlots/Anomalies/MozambiqueEddy19931995NotNorm_%d.png', i));
  writeVideo(writerObj, GeoThatDay);
end
close(writerObj);