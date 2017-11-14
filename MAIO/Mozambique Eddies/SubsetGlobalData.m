% Code for splitting the dataset so that we get a smaller subset that only
% contains the eddies within the Mozambique Channel

close all, clear all

% ncdisp('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc')

% import the data variables
% warning: this takes a long time so I suggest you run this as little as
% possible
fullData=zeros(22813018,9);
fullData(:,1)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','amplitude');
fullData(:,2)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','cyclonic_type');
fullData(:,3)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','latitude');
fullData(:,4)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','longitude');
fullData(:,5)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','observation_number');
fullData(:,6)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','speed_average');
fullData(:,7)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','speed_radius');
fullData(:,8)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','time');
fullData(:,9)=ncread('/Users/victoronink/Desktop/DataFiles/MAIO/eddy_trajectory_19930101_20160925.nc','track');

%%
%Define the domain of the Mozambique Channel
yv=[-10,-10,-12,-25.6,-26.5,-10];   %Degrees latitude south
xv=[32,40, 49.2, 45, 32,32];    %Degrees longitude East

%Define the new Mozambique Dataset
mozambiqueData=[];

N=1;
for index=0:max(fullData(:,9))
    if N+2200>length(fullData(:,9))
        k=find(index==fullData(N:end,9));
    else 
        k=find(index==fullData(N:N+2200,9));
    end
    inRange=inpolygon(fullData(N,4),fullData(N,3),xv,yv);
    if inRange==true
        mozambiqueData(end+1:end+length(k),:)=fullData(N+k-1,:);
    end
    N=N+length(k);
    index
end
%%
mozambiqueDataBackup=mozambiqueData;
%%
%save('MozamEddyDataTotal','mozambiqueData')