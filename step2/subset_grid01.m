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
geo_name='geo_em.d01_T3.nc';
% The ratio of WRF file resolution (93.75km here) to large-scale file
% resolution (250m here). So, 93.75km/250m = 375
div=375;
% First time running this code? For matching grids between the large-scale
% file and the WRF file
ifFirst = true;
%% read the large-scale file
input=ncreadall(input_name);
input_hgt=input.HGT_M;
input_hgt(isnan(input_hgt))=0;
load(inputLatLon);
% reshape the fine grid to 1 dimensional
% X=[reshape(input_lat,[],1),reshape(input_lon,[],1)];
%% read the WRF file
clear sub_hgt sub_lat sub_lon I_ind J_ind X
geo01=ncreadall(geo_name);
geo01_hgt=geo01.HGT_M;
geo01_lat=geo01.XLAT_M;
geo01_lon=geo01.XLONG_M;
% geo01_hgt(:,:)=0;
%% get the location for your subset
% % use the top corner grid as reference, find the closest grid to this
% reference in your large-scale file
if (ifFirst == true)
XI=[input_lat(1,1),input_lon(1,1)];
X=[reshape(geo01_lat,[],1),reshape(geo01_lon,[],1)];
% fast search of the location in the one dimensional file
k = dsearchn(X,XI);
% convert the location to 2 dimensional
[I1,J1]=ind2sub(size(geo01_lat),k);
XI=[input_lat(end,end),input_lon(end,end)];
k = dsearchn(X,XI);
[I2,J2]=ind2sub(size(geo01_lat),k);
save('matchPoints_dom1.mat','I1','J1','I2','J2');
else
load('matchPoints_dom1.mat');
end
[m, n]=size(geo01_lat);
[inm, inn]=size(input_lat);
% get the subset and put them into the three variables

for i=I1:I2
    for j=J1:J2
        if i~=I2 && j~=J2
I_ind=1+div*(i-I1):div*(i-I1+1);
J_ind=1+div*(j-J1):div*(j-J1+1);
        else
            I_ind=1+div*(i-I1):inm;
            J_ind=1+div*(j-J1):inn;
        end
sub_lat(i,j)=nanmean(nanmean(input_lat(I_ind,J_ind)));
sub_lon(i,j)=nanmean(nanmean(input_lon(I_ind,J_ind)));
sub_hgt(i,j)=nanmean(nanmean(input_hgt(I_ind,J_ind)));
    end
end
%% check if the latitudes do match well
sub_hgt(I2:m,:)=0;
contourf(sub_lat); % latitude cut from large-scale file
figure;contourf(geo01_lat); % latitude for the WRF file
figure;contourf(sub_hgt);  % result topography interpolated
%% replace the original height with the new one
ncid=netcdf.open(geo_name,'write');
ncwrite(geo_name,'HGT_M',sub_hgt);