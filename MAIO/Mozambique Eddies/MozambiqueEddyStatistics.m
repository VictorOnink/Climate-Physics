% Code for the computing of different characteristics of cyclonic and
% anticyclonic cyclones formed in all the different regions

close all
clear all

load('AllTracks.mat')

%% Divide all the tracks into whether they are cyclonic or anti cyclonic
NarrowCyc={};
NarrowAntiCyc={};
NorthCyc={};
NorthAntiCyc={};
SouthCyc={};
SouthAntiCyc={};

% Sort all the eddies from the narrows into cyclonic or anticyclonic
Cyc=1;
Acyc=1;
for i=1:length(BasinNarrow)
    track=BasinNarrow{i};
    if track(1,2)==1
        NarrowCyc{Cyc}=track;
        Cyc=Cyc+1;
    elseif track(1,2)==-1
        NarrowAntiCyc{Acyc}=track;
        Acyc=Acyc+1;
    end
end

% Sort all the eddies from the north basin into cyclonic or anticyclonic
Cyc=1;
Acyc=1;
for i=1:length(BasinNorth)
    track=BasinNorth{i};
    if track(1,2)==1
        NorthCyc{Cyc}=track;
        Cyc=Cyc+1;
    elseif track(1,2)==-1
        NorthAntiCyc{Acyc}=track;
        Acyc=Acyc+1;
    end
end

% Sort all the eddies from the south basin into cyclonic or anticyclonic
Cyc=1;
Acyc=1;
for i=1:length(BasinSouth)
    track=BasinSouth{i};
    if track(1,2)==1
        SouthCyc{Cyc}=track;
        Cyc=Cyc+1;
    elseif track(1,2)==-1
        SouthAntiCyc{Acyc}=track;
        Acyc=Acyc+1;
    end
end

figure(2)
for i=1:3
    %select which data we are dealing with first
   if i==1
       cyclonic=NorthCyc;
       anticyclonic=NorthAntiCyc;
   elseif i==2
       cyclonic=NarrowCyc;
       anticyclonic=NarrowAntiCyc;
   elseif i==3
       cyclonic=SouthCyc;
       anticyclonic=SouthAntiCyc;
   end
   
   %retrieve the cyclone data and for the amplitude put it into a long list
   %so that we can then plot the amplitude as a function of the eddy radius
   cyclone_amp=[];
   cyclone_time=[];
   for j=1:length(cyclonic)
       track=cyclonic{j};
       cyclone_amp(end+1:end+length(track(:,1)),:)=track(:,[1 7]);
       cyclone_time(end+1:end+length(track(:,1)),:)=track(:,[5 7]);
   end
   anticyclone_amp=[];
   anticyclone_time=[];
   for j=1:length(anticyclonic)
       track=anticyclonic{j};
       anticyclone_amp(end+1:end+length(track(:,1)),:)=track(:,[1 7]);
       anticyclone_time(end+1:end+length(track(:,1)),:)=track(:,[5 7]);
   end
   
   % Bin Averages
   dR=5; %km, bin size
   [BinMeanCycAmp, BinRadiusCycAmp]=BinAverage(dR, cyclone_amp(:,1), cyclone_amp(:,2) );
   [BinMeanACycAmp, BinRadiusACycAmp]=BinAverage(dR, anticyclone_amp(:,1), anticyclone_amp(:,2) );
   [BinMeanCycTime, BinRadiusCycTime]=BinAverage(dR, cyclone_time(:,1), cyclone_time(:,2) );
   [BinMeanACycTime, BinRadiusACycTime]=BinAverage(dR, anticyclone_time(:,1), anticyclone_time(:,2) );
   
   subplot(3,2,2*i-1)
       plot(BinRadiusCycAmp(~isnan(BinMeanCycAmp)),BinMeanCycAmp(~isnan(BinMeanCycAmp)),'r'), hold on
       plot(BinRadiusACycAmp(~isnan(BinMeanACycAmp)),BinMeanACycAmp(~isnan(BinMeanACycAmp)),'b')
   subplot(3,2,2*i)
       plot(BinRadiusCycTime(~isnan(BinMeanCycTime)),BinMeanCycTime(~isnan(BinMeanCycTime)),'r'), hold on
       plot(BinRadiusACycTime(~isnan(BinMeanACycTime)),BinMeanACycTime(~isnan(BinMeanACycTime)),'b')
end
subplot(3,2,1)
    title('Northern Basin')
    xlabel('radius (km)')
    ylabel('amplitude (cm)')
subplot(3,2,2)
    title('Northern Basin')
    xlabel('radius (km)')
    ylabel('life time (days)')
subplot(3,2,3)
    title('Channel Narrows')
    xlabel('radius (km)')
    ylabel('amplitude (cm)')
subplot(3,2,4)
    title('Channel Narrows')
    xlabel('radius (km)')
    ylabel('life time (days)')
subplot(3,2,5)
    title('Southern Basin')
    xlabel('radius (km)')
    ylabel('amplitude (cm)')
subplot(3,2,6)
    title('Southern Basin')
    xlabel('radius (km)')
    ylabel('life time (days)')