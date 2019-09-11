%% take wrf_hydro output, collect only gridded discharge files, rename them, rotate and clip them
close all; 
clear all;
%cd('~/HYDRO_analysis/matlab');

%% You ONLY need to edit this section to specify input/output filename, rotation angle and how you want to cut the domain after rotation
% specify the filename WRF-hydro 6hr-ly output, and output file name.
input_filename='201105241200.CHRTOUT_GRID3';
output_filenm=[input_filename,'_rot.nc'];
% rotation angle
angle=0;
% cut the file
left=257;
right=left+799;
top=101;
bottom=1220;
%%
display(input_filename);
tform = affine2d([cosd(angle) sind(angle) 0; -sind(angle) cosd(angle) 0; 0 0 1]);
allvar=ncreadall(input_filename); % read discharge variable from file
rows = right - left + 1;    
cols = bottom - top + 1;    
if exist(output_filenm,'file') == 2 
    eval(['delete ', output_filenm]);
end    
var = allvar.streamflow;    
ntime = size(var,3);    
RA=imref2d(size(var));   
% rotation
[var_rot, RB] = imwarp(var,RA,tform,'nearest');    
var_rot(var_rot<=0)= 1e33;    
% cut the domain out
var_clip=var_rot(left:right,top:bottom,:);    
%% write output file  
ncid = netcdf.create(output_filenm,'NETCDF4');    
dimcol = netcdf.defDim(ncid,'cols',rows);    
dimrow = netcdf.defDim(ncid,'rows',cols);    
dimt = netcdf.defDim(ncid,'Time',ntime);  
varid = netcdf.defVar(ncid,'streamflow','double',[dimcol dimrow dimt]);    
netcdf.defVarFill(ncid,varid,false,1e33);    
netcdf.close(ncid);    
ncwrite(output_filenm,'streamflow',var_clip);    
    
    
    
    
    




