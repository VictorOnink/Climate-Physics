% Code for the plotting of the cyclone tracks in different regions of the
% Mozambique Channel

close all
clear all

%Import the dataset, the variables inside the dataset are saved in the
%following columns 
%amplitude=mozambiqueData(:,1);
%type=mozambiqueData(:,2);
%lat=mozambiqueData(:,3);
%lon=mozambiqueData(:,4);
%obs_num=mozambiqueData(:,5);
%v_ave=mozambiqueData(:,6);
%v_r=mozambiqueData(:,7);
%t=mozambiqueData(:,8);
%track=mozambiqueData(:,9);
load('MozamEddyDataTotal.mat')

%Define a track
Tracks={};
k=1;    %track number
j=1;    %element number
for i=1:length(mozambiqueData(:,9))
    if i==1
        iden=mozambiqueData(i,9);
    else
        if mozambiqueData(i,9)==iden
            continue
        else
            datatrack=mozambiqueData(j:j+(i-j-1),:);
            Tracks{k}=datatrack;
            k=k+1;
            j=i;
            iden=mozambiqueData(i,9);
        end
    end
end
%% Basin Definitions for the sorting of the starting points of the eddies
%Define the total domain of the Mozambique Channel
yv_tot=[-10,-10,-12,-25.6,-26.5,-10];   %Degrees latitude south
xv_tot=[32,40, 49.2, 45, 32,32];    %Degrees longitude East
%Northern basin of the Mozambique Channel
yv_nor=[-10,-12,-16.5,-14.5,-10];
xv_nor=[40,49.2,46,40,40];
%Narrows of the Mozambique Channel
yv_nar=[-14.5,-16.5,-18.5,-17,-14.5];
xv_nar=[40,46,44.3,38.5,40];
%Southern basin of the Mozambique Channel
yv_sou=[-17,-18.5,-25.6,-26.5,-20,-17];
xv_sou=[38.5,44.3,45,32,32,38.5];

%% Break the tracks into which basin they start out in
BasinNorth={};
BasinNarrow={};
BasinSouth={};

Nor=1;
Nar=1;
Sou=1;

for index=1:length(Tracks)
    track=Tracks{index};
    RangeNorth=inpolygon(track(1,4),track(1,3),xv_nor,yv_nor);
    if RangeNorth==true      
        BasinNorth{Nor}=track;
        Nor=Nor+1;     
    else
        RangeSouth=inpolygon(track(1,4),track(1,3),xv_sou,yv_sou);
        if RangeSouth==true           
            BasinSouth{Sou}=track;
            Sou=Sou+1;
        else
            RangeNarrow=inpolygon(track(1,4),track(1,3),xv_nar,yv_nar);
            if RangeNarrow==true           
                BasinNarrow{Nar}=track;
                Nar=Nar+1;             
            end
        end
    end
    disp(index)
end

save('AllTracks','Tracks','BasinNorth','BasinNarrow','BasinSouth')
%% Break separated tracks by lifetime
agesplit=60; %days
agesplit2=120;
NorthYoung=[];
NorthOld=[];
NorthOldOld=[];
young=1;
old=1;
oldold=1;

for i=1:length(BasinNorth)
    track=BasinNorth{i};
    if length(track)<=agesplit
        NorthYoung{young}=track;
        young=young+1;
    elseif length(track)<=agesplit2
        NorthOld{old}=track;
        old=old+1;
    else
        NorthOldOld{oldold}=track;
        oldold=oldold+1;
    end
end

NarrowYoung=[];
NarrowOld=[];
NarrowOldOld=[];
young=1;
old=1;
oldold=1;

for i=1:length(BasinNarrow)
    track=BasinNarrow{i};
    if length(track)<=agesplit
        NarrowYoung{young}=track;
        young=young+1;
    elseif length(track)<=agesplit2
        NarrowOld{old}=track;
        old=old+1;
    else
        NarrowOldOld{oldold}=track;
        oldold=oldold+1;
    end
end

SouthYoung=[];
SouthOld=[];
SouthOldOld=[];
young=1;
old=1;
oldold=1;

for i=1:length(BasinSouth)
    track=BasinSouth{i};
    if length(track)<=agesplit
        SouthYoung{young}=track;
        young=young+1;
    elseif length(track)<=agesplit2
        SouthOld{old}=track;
        old=old+1;
     else
        SouthOldOld{oldold}=track;
        oldold=oldold+1;    
    end
end

%% plotting on a map of the Mozambique channel subplots
% split by anti/cyc
% 1 = red = anticyclonic
% split by age
for i=1:3
    figure(i)
    subplot(2,3,1)
    title(['Anticyclones <',num2str(agesplit),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','off','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','on','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');    
    if i==1
        subtracks=NorthYoung;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowYoung;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthYoung;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','red',...
                'DisplayName','Anti-Cyclonic');
            hold on
        end
    end
end
for i=1:3
    figure(i)
    subplot(2,3,4)
    title(['Cyclones <',num2str(agesplit),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','on','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','on','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');
    if i==1
        subtracks=NorthYoung;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowYoung;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthYoung;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==-1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','blue',...
                'DisplayName','Cyclonic');
            hold on
        end
    end
end
for i=1:3
    figure(i)
    subplot(2,3,2)
    title(['Anticyclones >',num2str(agesplit),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','off','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','off','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');
    if i==1
        subtracks=NorthOld;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowOld;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthOld;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','red',...
                'DisplayName','Anti-Cyclonic');
            hold on
        end
    end
end
for i=1:3
    figure(i)
    subplot(2,3,5)
    title(['Cyclones >',num2str(agesplit),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','on','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','off','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');
    if i==1
        subtracks=NorthOld;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowOld;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthOld;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==-1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','blue',...
                'DisplayName','Cyclonic');
            hold on
        end
    end
end
for i=1:3
    figure(i)
    subplot(2,3,3)
    title(['Anticyclones >',num2str(agesplit2),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','off','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','off','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');
    if i==1
        subtracks=NorthOldOld;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowOldOld;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthOldOld;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','red',...
                'DisplayName','Anticyclonic');
            hold on
        end
    end
end
for i=1:3
    figure(i)
    subplot(2,3,6)
    title(['Cyclones >',num2str(agesplit2),' days'])
    map=axesm('miller','maplatlim',[-40, -5],'maplonlim',[15, 52],'grid','on',...
        'meridianlabel','on','mlabellocation',[20 25 30 35 40 45 50],...
        'mlabelparallel','south',...
        'parallellabel','off','plabellocation',[-10 -15 -20 -25 -30 -35]);    
    setm(map,'mlinelocation',5,'plinelocation',5)
    tightmap; 
    load coastlines
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');
    if i==1
        subtracks=NorthOldOld;
        plotm(yv_nor,xv_nor,'LineWidth',1,'Color','black')
        hold on
    elseif i==2
        subtracks=NarrowOldOld;
        plotm(yv_nar,xv_nar,'LineWidth',1,'Color','black')
        hold on
    elseif i==3
        subtracks=SouthOldOld;
        plotm(yv_sou,xv_sou,'LineWidth',1,'Color','black')
        hold on
    end
    for i=1:length(subtracks)
        track=subtracks{i};
        if track(1,2)==-1
            plotm(track(:,3),track(:,4),'LineWidth',1,'Color','blue',...
                'DisplayName','Cyclonic');
            hold on
        end
    end
end
%% Birthplace. Actually useful?
figure()
title('Eddy Birthplace')
map=axesm('miller','maplatlim',[-30, -10],'maplonlim',[30, 52],'grid','on',...
    'meridianlabel','on','mlabellocation',[ 35 40 45 50],...
    'mlabelparallel','south',...
    'parallellabel','on','plabellocation',[-15 -20 -25]);    
setm(map,'mlinelocation',5,'plinelocation',5)
tightmap; 
load coastlines
geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','Green');

for i=1:length(mozambiqueData)
    if mozambiqueData(i,5)==1 % First point. change this to 0 later
        if mozambiqueData(i,2)==1 % Anti
            plotm(mozambiqueData(i,3),mozambiqueData(i,4),'rx')
            hold on
        else % Cyclone
            plotm(mozambiqueData(i,3),mozambiqueData(i,4),'bx')
            hold on
        end
    else
    end
end       
%% Save the Tracks of the different basins
% save('AllTracks','Tracks','BasinNorthCyc','BasinNarrowCyc','BasinSouthCyc','BasinNorthAnti','BasinNarrowAnti','BasinSouthAnti')
