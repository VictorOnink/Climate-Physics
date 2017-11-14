function [curlz,cav] = curlm(lat,lon,U,V)
%CURLM computes the curl of U and V given their geographical coordinates
%lat and lon.  lat, lon, U, and V must be m by n matrices of the same size.
% Coordinates given by lat and lon must be regular grids of latitudes and
% longitudes--that is, they must be separated as 120E, 121E, 123E, etc, not
% equally spaced in terms of kilometers of separation.  
% 
% Requires mapping toolbox. 
% Chad Greene, Jan 2014.

latsize = size(lat); 



x_m = zeros(latsize); % preallocate x_m
dx_m = distance('rh',lat(:,1),lon(:,2),lat(:,1),lon(:,1),6378137); % x distance between grid points
dx_m_grid = repmat(dx_m,1,latsize(2)-1); 

x_m(:,2:end) = cumsum(dx_m_grid,2); % x coordinates in meters, origin to the left


y_m = zeros(latsize); 
dy_m = distance('rh',lat(2,:),lon(1,:),lat(1,:),lon(1,:),6378137);
dy_m_grid = repmat(dy_m,latsize(1)-1,1); 
y_m(2:end,:) = cumsum(dy_m_grid,1); 

if lat(2,1)<lat(1,1)
    y_m = flipud(y_m); 
end

[curlz,cav] = curl(x_m,y_m,U,V); 


end