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
geo_name='geo_em.d02_T3.nc';
% The ratio of WRF file resolution (18.75km here) to large-scale file
% resolution (250m here). So, 18.75km/250m = 375
div=75;
% First time running this code? For matching grids between the large-scale
% file and the WRF file
ifFirst = true;
%%  read the large-scale file
load(inputLatLon);
input=ncreadall(input_name);
input_hgt=input.HGT_M;
input_hgt(isnan(input_hgt))=0;
%% read the WRF file
clear sub_hgt sub_lat sub_lon I_ind J_ind
geo02=ncreadall(geo_name);
geo02_hgt=geo02.HGT_M;
geo02_lat=geo02.XLAT_M;
geo02_lon=geo02.XLONG_M;
%% get the location and sebset
if (ifFirst == true)
% % use the top corner grid as reference
X=[reshape(input_lat,[],1),reshape(input_lon,[],1)];
XI=[geo02_lat(1,1),geo02_lon(1,1)];
%  fast search of the location in the one dimensional file
k = dsearchn(X,XI);
%  convert the location to 2 dimensional
[I,J]=ind2sub(size(input_lat),k);
save('matchPoints_dom2.mat','I','J');
else
% I=3038;
% J=28538;
load('matchPoints_dom3.mat');
end
% find the size of the small subset
[m, n]=size(geo02_lat);
% get the subset and put them into the three variables
I=I-div/2;
J=J-div/2;
for i=1:m
    for j=1:n
I_ind=I+div*(i-1):I+div*i-1;
J_ind=J+div*(j-1):J+div*j-1;
sub_lat(i,j)=nanmean(nanmean(input_lat(I_ind,J_ind)));
sub_lon(i,j)=nanmean(nanmean(input_lon(I_ind,J_ind)));
sub_hgt(i,j)=nanmean(nanmean(input_hgt(I_ind,J_ind)));
    end
end
%% check if the latitudes do match well
contourf(sub_lat);  % latitude cut from large-scale file
figure;contourf(geo02_lat); % latitude for the WRF file
figure;contourf(sub_hgt);  % result topography interpolated
%% replace the original height with the new one
ncid=netcdf.open(geo_name,'write');
ncwrite(geo_name,'HGT_M',sub_hgt);