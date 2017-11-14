%  Halo Fig 8, Density Distribution as function of Radius
% Take the average radius of the eddy as qualifier
% ...and other stats bits Halo fig 10
close all
clear all

load('AllTracks.mat')

% For checking ecsapees
%Define the total domain of the Mozambique Channel
yv_tot=[-10,-10,-12,-25.6,-26.5,-10];   %Degrees latitude south
xv_tot=[32,40, 49.2, 45, 32,32];    %Degrees longitude East

for j=1:4
    % Select data to work on
    if j==1
        stuff=Tracks;
    elseif j==2
        stuff=BasinNarrow;
    elseif j==3
        stuff=BasinNorth;
    elseif j==4
        stuff=BasinSouth;
    end
    
    % Empties to fill
    meanradtot={};
    meanamptot={};
    lifetot={};
    disttot={};
    meanradAnti={};
    meanampAnti={};
    lifeAnti={};
    distAnti={};
    meanradCyc={};
    meanampCyc={};
    lifeCyc={};
    distCyc={};
    tot=1;
    anti=1;
    cyc=1;
    
    % Calculate, for each eddy in each section
    for i=1:length(stuff)
        eddy=stuff{i};
        [l,w]=size(eddy);
        if l~=1 % Filtering out singles. Will be obselete when we re-run selection        
            meanradtot{tot}=mean(eddy(:,7));
            meanamptot{tot}=mean(eddy(:,1));
            lifetot{tot}=length(eddy(:,8));
            travel=0;
            for z=1:length(eddy)-1
                [Haversine,~]=lldistkm([eddy(z,3) eddy(z,4)],[eddy(z+1,3) eddy(z+1,4)]);
                travel=travel+Haversine;
            end
            
            disttot{tot}=travel;
            tot=tot+1;
            if eddy(:,2)==1 % Antis
                meanradAnti{anti}=mean(eddy(:,7));
                meanampAnti{anti}=mean(eddy(:,1));
                lifeAnti{anti}=length(eddy(:,8));
                travel=0;
                for z=1:length(eddy)-1
                    [Haversine,~]=lldistkm([eddy(z,3) eddy(z,4)],[eddy(z+1,3) eddy(z+1,4)]);
                    travel=travel+Haversine;
                end  
                distAnti{anti}=travel;
                anti=anti+1;
            else % Cyclones
                meanradCyc{cyc}=mean(eddy(:,7));
                meanampCyc{cyc}=mean(eddy(:,1));
                lifeCyc{cyc}=length(eddy(:,8));
                travel=0;
                for z=1:length(eddy)-1
                    [Haversine,~]=lldistkm([eddy(z,3) eddy(z,4)],[eddy(z+1,3) eddy(z+1,4)]);
                    travel=travel+Haversine;
                end
                distCyc{cyc}=travel;
                cyc=cyc+1;
            end
        end
    end
    % Convert cells to matrix format for plottability
    meanradtot=cell2mat(meanradtot);
    meanamptot=cell2mat(meanamptot);
    lifetot=cell2mat(lifetot);
    disttot=cell2mat(disttot);
    meanradAnti=cell2mat(meanradAnti);
    meanampAnti=cell2mat(meanampAnti);
    lifeAnti=cell2mat(lifeAnti);
    distAnti=cell2mat(distAnti);
    meanradCyc=cell2mat(meanradCyc);
    meanampCyc=cell2mat(meanampCyc);
    lifeCyc=cell2mat(lifeCyc);
    distCyc=cell2mat(distCyc);

% And potting all the things. Literally. All of them. How many subplots is
% too many subplots? Keep travels separate for now
figure(j)
screensize = get( groot, 'Screensize' );
set(figure(j),'Position',screensize)
% Density distributions
subplot(3,3,1)
    tothist=histogram(meanradtot,'BinWidth',10,'Normalization','probability','Facecolor','black');
    if j==1
        title('Total')
    elseif j==2
        title('Narrow')
    elseif j==3
        title('North')
    else
        title('South')
    end
    ylabel('Eddies (%)')
    if j==3
        Denseax=gca; % Because this one has the biggest range. Put everything on the same
    end
subplot(3,3,2)
    antichrist=histogram(meanradAnti,'NumBins',tothist.NumBins,'BinWidth',tothist.BinWidth);
    histogram('BinEdges',antichrist.BinEdges,'BinCounts',antichrist.Values/length(tothist.Data),'Facecolor','red')
    title('Anticyclones')
    yticklabels([])
    grid on
    grid minor
subplot(3,3,3)
    cychist=histogram(meanradCyc,'NumBins',tothist.NumBins,'BinWidth',tothist.BinWidth);
    histogram('BinEdges',cychist.BinEdges,'BinCounts',cychist.Values/length(tothist.Data),'Facecolor','blue')
    title('Cyclones')
    yticklabels([])
% Amplitudes
subplot(3,3,4)
    [Ampmean, Binradius]=BinAverage(10,meanamptot,meanradtot);
    plot(Binradius,Ampmean,'k-','LineWidth',2)
    hold on
    plot(meanradtot,meanamptot,'k.')
    if j==1
        ampax=gca;
    end
    ylim(ampax.YLim)
    ylabel('Amplitude (cm)')
subplot(3,3,5)
    [Ampmean, Binradius]=BinAverage(10,meanampAnti,meanradAnti);
    plot(Binradius,Ampmean,'r-','LineWidth',2)
    hold on
    plot(meanradAnti,meanampAnti,'r.')
    ylim(ampax.YLim)
    yticklabels([])
subplot(3,3,6)
    [Ampmean, Binradius]=BinAverage(10,meanampCyc,meanradCyc);
    plot(Binradius,Ampmean,'b-','LineWidth',2)
    hold on
    plot(meanradCyc,meanampCyc,'b.')
    ylim(ampax.YLim)
    yticklabels([])

% Lifetimes
subplot(3,3,7)
    [Lifemean, Binradius]=BinAverage(10,lifetot,meanradtot);
    plot(Binradius,Lifemean,'k-','LineWidth',2)
    hold on
    plot(meanradtot,lifetot,'k.');
    ylabel('Lifetime (days)')
subplot(3,3,8)
    [Lifemean, Binradius]=BinAverage(10,lifeAnti,meanradAnti);
    plot(Binradius,Lifemean,'r-','LineWidth',2)
    hold on
    plot(meanradAnti,lifeAnti,'r.')
    yticklabels([])
    xlabel('Radius (km)')
subplot(3,3,9)
    [Lifemean, Binradius]=BinAverage(10,lifeCyc,meanradCyc);
    plot(Binradius,Lifemean,'b-','LineWidth',2)
    hold on
    plot(meanradCyc,lifeCyc,'b.');
    yticklabels([])
    
% Travel distances
figure(j+4)
subplot(3,2,1)
    plot(meanradtot,disttot,'k.')
    if j==1
        title('Total')
    elseif j==2
        title('Narrow')
    elseif j==3
        title('North')
    else
        title('South')
    end
    xlabel('Mean radius (km)')
    ylabel('Travel Distance (km)')
subplot(3,2,2)
    plot(lifetot,disttot,'k.')
    xlabel('Lifetime (days)')
subplot(3,2,3)
    plot(meanradAnti,distAnti,'r.')
subplot(3,2,4)
    plot(lifeAnti,distAnti,'r.')
subplot(3,2,5)
    plot(meanradCyc,distCyc,'b.')
subplot(3,2,6)
    plot(lifeCyc,distCyc,'b.')
    
    
end
% Evening out the axes 
for j=1:4
    figure(j)
    for i=1:3
        subplot(3,3,i)
        yticks(Denseax.YTick)
        if i==1
            yticklabels(Denseax.YTick*100)
        end
        ylim(Denseax.YLim)
    end
    for i=1:9
        subplot(3,3,i)
        xlim([50 200])
        grid on
        grid minor
    end
    for i=7:9
        subplot(3,3,i)
        ylim([0 500])
    end
end
