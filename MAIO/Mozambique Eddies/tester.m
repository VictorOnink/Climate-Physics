% Finding which Narrow eddies exist purely within the specified period, by
% which we can match to geostrophic data
clear all
close all

% load('AllTracks.mat')

% % Data period already downloaded
% start=15706;
% stop=16436;
% match={};
% m=1;
% 
% for i=1:length(BasinNarrow)
%     eddy=BasinNarrow{i};
%     if eddy(1,8)>=start && eddy(end,8)<=stop
%         match{m}=eddy;
%         m=m+1;
%     end      
% end

% % Only have three narrow eddies in that period. Different geostrophic
% % dataset?
% % If we restrict ourselves to downloading a 2-year set, best coverage?
% year={};
% for i=1:length(BasinNarrow)
%     eddy=BasinNarrow{i};
%     if length(eddy)<=730
%         year{i}=[eddy(1,8),eddy(end,8)];
%     end
% end

% % Ok they all less than 2 years long
% % Optimise to find most in a 2-year set
% saver={};
% len=0;
% jj=1;
% while jj<length(BasinNarrow)-len
%     starter=BasinNarrow{jj};
%     start=starter(1,8);
%     stop=start+730;
%     matches={};
%     mat=1;
%     for i=jj:length(BasinNarrow)
%         eddy=BasinNarrow{i};
%         if eddy(end,8)<stop
%             matches{mat}=eddy(:,:);
%             mat=mat+1;
%         else
%             break
%         end
%     end
%     if length(matches)>len
%         len=length(matches);
%         saver=matches;
%         beststart=start;
%         beststop=stop;
%     end
%     jj=jj+1;
% end
            
%% Better Bathymetric mapping
% START CODE
addpath(genpath('/Users/victoronink/Documents/MATLAB/m_map'))


% define below the coordinates of your framework
    lon_min = 15; lon_max = 52;
    lat_min  = -40;  lat_max = -5;

% extract bathymetry for the region of study    
    [elevation, vlat, vlon] = mygrid_sand_1m([lat_min lat_max lon_min lon_max],1);
%     vlon    = vlon - 360; % you may need to comment this line depending on your longitudes of study
    
% initiate a map projection
    m_proj('mercator','lat',[lat_min lat_max],'lon',[lon_min lon_max]); hold on
    
% plot bathymetry
    [c1,h1] = m_contourf(vlon,vlat,elevation,[-1000 -2000 -3000 -4000 -5000 -6000]); hold on
    set(h1,'color','k');
    colormap('gray')
    c=colorbar;
    c.Label.String = 'Depth contours (m)';

    
% save full resolution coastline under desired path and filename (define below)
%    m_gshhs_f('save','/Users/victoronink/Desktop/DataFiles/MAIO/moz_coast');
% comment the above line after run it one time since your coastline is then created and saved. You don't need to have it active everytime, it requires lot (several seconds) of computing time.

% call the coastline previously saved
    m_usercoast('moz_coast.mat','patch',[.7 .7 .7]);
    m_grid('linestyle','none','tickedir','out','fontsize',14);


% END CODE
        
        
        
        
        
        
        
        


        