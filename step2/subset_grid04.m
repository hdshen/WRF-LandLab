%%
clear;
close all;
%% specify your input file names. You ONLY need to change this part to process your own files
% This is the large-scale file with the entire mountain range
input_name= '../geo_em.coor_gradient_100m_T3.nc';
% This file contains the coordinates for the large-scale file, 
% with the coordinates name as input_lat and input_lon, and each of two varirale 
% has the same dimension as the large-scale file
inputLatLon='../latlon.mat';
% This is the file name for your nested domain WRF file
geo_name='geo_em.d04_T3.nc';
% The ratio of WRF file resolution (250 m here) to large-scale file
% resolution (250m here). So, 250m250m = 1
div=1;
% First time running this code? For matching grids between the large-scale
% file and the WRF file
ifFirst = true;
%%  read the large-scale file
input=ncreadall(input_name);
input_hgt=input.HGT_M;
input_hgt(isnan(input_hgt))=0;

%% read the WRF file
clear sub_hgt sub_lat sub_lon I_ind J_ind
geo04=ncreadall(geo_name);
geo04_hgt=geo04.HGT_M;
geo04_lat=geo04.XLAT_M;
geo04_lon=geo04.XLONG_M;
%% get the location and sebset
if (ifFirst == true)
% reshape the fine grid to 1 dimensional
X=[reshape(input_lat,[],1),reshape(input_lon,[],1)];
% use the top corner grid as reference
XI=[geo04_lat(1,1),geo04_lon(1,1)];
% fast search of the location in the one dimensional file
k = dsearchn(X,XI);
% convert the location to 2 dimensional
[I,J]=ind2sub(size(input_lat),k);
save('matchPoints_dom4.mat','I','J');
else
% I=5251;
% J=30151;
load('matchPoints_dom4.mat');
end
% find the size of the small subset
[m, n]=size(geo04_lat);
% get the subset and put them into the three variables
I_ind=I:I+m-1;
J_ind=J:J+n-1;
sub_lat=input_lat(I_ind,J_ind);
sub_lon=input_lon(I_ind,J_ind);
sub_hgt=input_hgt(I_ind,J_ind);
%% check if the latitudes do match well
contourf(sub_lat); % latitude cut from large-scale file
figure;contourf(geo04_lat); % latitude for the WRF file
figure;contourf(sub_hgt);   % result topography interpolated
%% replace the original height with the new one
ncid=netcdf.open(geo_name,'write');
ncwrite(geo_name,'HGT_M',sub_hgt);