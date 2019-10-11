% Written by Brigid Lynch (Indiana University, Dept. of Earth and
% Atmospheric Sciences. This code is a shortened example to show how to rotate 
% two Landlab topographies to the angle of the coastline, place them into the 
% large-scale topography file and interpolate between the two to build a 
% mountain range. To build a larger mountain range, you would repeat this 
% process with any number of Landlab Domains.

% NOTE: This script uses the gIDW function, included in this example folder and
% available on Mathworks File Exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/27562-inverse-distance-weighted-idw-or-simple-moving-average-sma-interpolation?focused=5152732&tab=function

close all
clear all
%% Rotate Landlab topographies %%

% We first need the individual Landlab topographies rotated to the desired
% angle to fit along the coast line. For this example: input_topo_1 needs 
% no rotation and input_topo_2 needs to be rotated -27 deg. The rotation 
% angles were determined manually for each domain based on their location 
% along the coastline.

% Read in the topography created in Landlab from the driver_create_topo.py
% file.

input_topo_n=ncread('landlab_topo_n.nc', 'topographic__elevation'); 

% All of the edges of the topography created in Landlab go to ~zero, but we
% want the North and South Edges (800 cells) to have topography that
% matches the mountain range, so we copy the neighboring elevations to
% achieve this. 

input_topo_n(:,1) = input_topo_n(:,2);
input_topo_n(:,end) = input_topo_n(:,end-1);

%%% Landlab Topography 1 %%%

% Since input_topo_1 is not rotated, we can keep it in it's original
% orientation.
input_topo_1 = input_topo_n;

%%% Landlab Topography 2 %%%
% This topography needs to be rotated to fit the angle of the coastline.

theta = -27; % Angle to rotate the topography     

% Perform rotation
rmat = [
cosd(theta) sind(theta) 0
-sind(theta) cosd(theta) 0
0           0          1];

mx = size(input_topo_n,2);
my = size(input_topo_n,1);
corners = [
    0  0  1
    mx 0  1
    0  my 1
    mx my 1];
new_c = corners*rmat;

T = maketform('affine', rmat);   

input_topo_2 = imtransform(input_topo_n, T, ...
    'XData',[min(new_c(:,1)) max(new_c(:,1))],...
    'YData',[min(new_c(:,2)) max(new_c(:,2))]);


input_topo_2 = input_topo_2(1:end-1,2:end); %trim edges with only zeros


%% ADD LAND/OCEAN TO LARGE-SCALE TOPOGRAPHY FILE %%
% Use a landmask file for the desired locations (Here: South America) to
% add land and ocean to the large-scale file (geo_em.coor.nc). 

% NOTE: We used the landmask function available on Matlab File Exchange to 
% build this file. 
% https://www.mathworks.com/matlabcentral/fileexchange/48661-landmask


% load landmask 
land=ncread('South_America_landmask.nc','landmask'); % 0=ocean; 1=land

% load elevation variable from large-scale file
topo=ncread('geo_em.coor.nc','HGT_M');

% set elevation of land to 0 and ocean to NaN in elevation variable
for i=1:size(land,1)
    for j=1:size(land,2)
        if land(i,j)==0
            topo(i,j)=NaN;
        elseif land(i,j)==1
            topo(i,j)=0;
        end
    end
end


%% Place rotated Landlab topographies onto large-scale topography file
% NOTE: this requires manual placement of each topography, for this project
% we placed the topographies along the Pacific coast of South America and
% the numbers used here, as well as the rotation of each topography, are
% specific to that location.

% 1st Landlab topography
for k=1:size(input_topo_1,1)
    for l=1:size(input_topo_1,2)
        % add the elevation variable from the first Landlab topography to
        % the large-scale topography variable (only add where there is
        % land, not ocean)
        if ~isnan(topo(k+5506,l+30250))
            topo(k+5506,l+30250)=input_topo_1(k,l);
        elseif isnan(topo(k+5506,l+30250))
            
        end
    end
end

% 2nd Landlab topography
for d=1:size(input_topo_2,1)
    for e=1:size(input_topo_2,2)
        if ~isnan(topo(d+6200,e+25700))
                    topo(d+6200,e+25700)=input_topo_2(d,e);
        elseif  isnan(topo(d+6200,e+25700))
                
        end
    end
end

%% Record the coordinates of the left and right side of each topo 
% (defined as the side with 800 cells)
% First do this within the individual topography variable (input_topo_1/2), 
% then translate this to the coordinates within large-scale topography file.
% For this example, we will be connected the left side of input_topo_1 to
% the right edge of input_topo_2

% NOTE: Coordinates within large-scale topography file must be selected
% manually as they will depend on the placement of individual Landlab
% topographies. 

% input_topo_1 is not rotated so this step is
% different for input_topo_1 than input_topo_2           

% LEFT SIDE OF TOPO 1
rows_1_L = 1:800;  % coordinates within input_topo_1 variable
for i=1:800
    % coordinates within large-scale topography file
    cols_1_L_TOPO(i)=30251;
    rows_1_L_TOPO(i)=rows_1_L(i)+5506;
end


% RIGHT SIDE OF TOPO 2
% coordinates within input_topo_2 variable
for j=1:size(input_topo_2,2) % count through every column of input input_topo_2
  R2 = find(input_topo_2(:,j), 1, 'first'); % for all rows, find the first 
                           % nonzero cell in each column and record its row
    if ~isempty(R2) 
      rows_2_R(j)=R2; % record all rows - here you are recording the row 
      % coorindate of both the top and right edges of the topography
      % in the next step you will limit this to only the right edge of the
      % topography
    end
    
end

for i=998:1362 % must manually find the columns of right edge of input_topo_2
    % translate coordinates within large-scale topography file
    cols_2_R_TOPO(i-997)=i+25700; 
    rows_2_R_TOPO(i-997)=rows_2_R(i)+6200;
    
end

%% Record the elevation values of the ridgeline of each topo - 
% separated into right and left sides of each topography

% Do this within the individual topography variable (input_topo_1/2)
% NOTE: coordinates that separate right and left ridgeline must be selected 
% manually as they will depend on the placement and rotation of individual 
% Landlab topographies.


for i=1:size(input_topo_1,2)-560 % run through half of the columns in input_topo_1
   
    ridgeline_1_L(i)=max(input_topo_1(:,i)); %record the maximum elevation of input_topo_1 across each column - left side
    avg_elev_1_L=mean(ridgeline_1_L); %find average elevation of the left ridegline
    ridgeline_1_R(i)=max(input_topo_1(:,i+560));%record the maximum elevation of input_topo_1 across each column - right side
    avg_elev_1_R=mean(ridgeline_1_R);%find average elevation of the right ridegline
end

for bb=1:size(input_topo_2,2)-866 % must manually select the length of the ridgeline depending on topography rotation
    ridgeline_2_L(bb)=max(input_topo_2(:,bb+184)); %record the maximum elevation of topo2 across each column - left side
    avg_elev_2_L=mean(ridgeline_2_L);%find average elevation of the left ridegline
    ridgeline_2_R(bb)=max(input_topo_2(:,bb+681));%record the maximum elevation of topo2 across each column - right side
    avg_elev_2_R=mean(ridgeline_2_R);%find average elevation of the right ridegline
end

%% Record the coordinates of the top and bottom side of each topo (defined as the side with 1120 cells)
% First do this within the individual topography (input_topo_1/2) then
% translate within the large-scale topography file

% NOTE: coordinates within large-scale topography file must be selected
% manually as they will depend on the placement of individual Landlab
% topographies. 

% Top & Bottom of input_topo_1
%coordinates within input_topo_1 variable 
for i=1:size(input_topo_1,2)
    T1_top = find(input_topo_1(:,i), 1, 'first');% for all rows, find the first nonzero cell in each column and record its row
    if ~isempty(T1_top)
      rows_top_1(i)=T1_top; %record this row coordinate
    end
    T1_bot = find(input_topo_1(:,i),1,'last');% for all rows, find the last nonzero cell in each column and record its row
    if ~isempty(T1_bot)
      rows_bot_1(i)=T1_bot;  %record this row coordinate
    end
end

for i=1:size(input_topo_1,2)
    side_1_top_TOPO(i)=input_topo_1(rows_top_1(i), i); % record the elevation values of the top of input_topo_1 in new variable
    
    %coordinates within large-scale topography
    rows_1_top_TOPO(i)=rows_top_1(i)+5506;
    cols_1_top_TOPO(i)=i+30250;
    
    avg_elev_1_top=mean(side_1_top_TOPO); % find the average elevation of top side of input_topo_1
    
    side_1_bot_TOPO(i)=input_topo_1(rows_bot_1(i), i); % record the elevation values of the bottom of input_topo_1 in new variable
    
    %coordinates within large-scale topography
    rows_1_bot_TOPO(i)=rows_bot_1(i)+5506;
    cols_1_bot_TOPO(i)=i+30250;
    
    avg_elev_1_bot=mean(side_1_bot_TOPO);% find the average elevation of bottom side of input_topo_1
end

% Top and bottom of input_topo_2
%coordinates within input_topo_2 variable 
for i=1:size(input_topo_2,2)
    T2_top = find(input_topo_2(:,i),1,'first'); %for all rows, find the first nonzero cell in each column and record its row
    if ~isempty(T2_top)
        rows_top_2(i)=T2_top;%record this row coordinate, here you are recording the row 
      % coorindate of both the top and right edges of the topography
      % in the next step you will limit this to only the top edge of the
      % topography
    end
    T2_bot = find(input_topo_2(:,i),1,'last');% for all rows, find the last nonzero cell in each column and record its row
    if ~isempty(T2_bot)
        rows_bot_2(i)=T2_bot;%record this row coordinate, here you are recording the row 
      % coorindate of both the bottom and left edges of the topography
      % in the next step you will limit this to only the bottom edge of the
      % topography
    end
end
   
for i=1:998
    side_2_top_TOPO(i)=input_topo_2(rows_top_2(i),i);% record the elevation values of the top of input_topo_2 in new variable
    %coordinates within large-scale topography
    rows_2_top_TOPO(i)=rows_top_2(i)+6200;
    cols_2_top_TOPO(i)=i+25700;
    
    avg_elev_2_top=mean(side_2_top_TOPO); % find the average elevation of top side of input_topo_2

    side_2_bot_TOPO(i)=input_topo_2(rows_bot_2(i+364),i+364);  % record the elevation values of the bottom of input_topo_2 in new variable
    
    %coordinates within large-scale topography
    rows_2_bot_TOPO(i)=rows_bot_2(i+364)+6200;
    cols_2_bot_TOPO(i)=i+26064;
    
    avg_elev_2_bot=mean(side_2_bot_TOPO);% find the average elevation of bottom side of input_topo_2
end

%% Adding average elevations to points between Landlab topographies, top and bottom edges
% In this step, we are first adding points of average elevations betweent
% the two topographies, then adding lines of elevation values that connected 
% the top and bottom edges of the two topographies. 

% Edge elevation half way between input_topo_1 and input_topo_2, top

% Find the mid point row/col for the top edge and insert average elevation
% within the large-scale topography file.
mid_point_1_top_row = round((rows_1_top_TOPO(1)+rows_2_top_TOPO(end))/2); %find row
mid_point_1_top_col = round((cols_1_top_TOPO(1)+cols_2_top_TOPO(end))/2); %find column
mid_point_1_top_elev = (avg_elev_1_top+avg_elev_2_top)/2; %find elevation
topo(mid_point_1_top_row,mid_point_1_top_col)=mid_point_1_top_elev; %place within large-scale topography variable

% Create edge line  between input_topo_1 and midpoint, top

% Create variables with rows/cols of top of input_topo_1 and elevation at these
% locations
Xc_1_mid_top = [rows_1_top_TOPO(1),mid_point_1_top_row];% top row of input_topo_1 and row of mid point
Yc_1_mid_top = [cols_1_top_TOPO(1),mid_point_1_top_col ];% top column of input_topo_1 and column of mid point
Vc_1_mid_top = [avg_elev_1_top, mid_point_1_top_elev];% average elevation of top of input_topo_1 and mid point

% Create a line with 800 points using the coordinates from above.
[fill_rows_1_mid_top, fill_cols_1_mid_top, I_1_mid_top]=improfile(topo,Xc_1_mid_top,Yc_1_mid_top,800);

% Use the coordinates from this line along with the coordinates and
% elevations of the top endpoints of input_topo_1 and midpoint to calculate
% elevations along this line. We use the gIDW which computes the inverse
% weighted interpolation. 
coords_1_mid_top=round([fill_rows_1_mid_top, fill_cols_1_mid_top]);
Xi_1_mid_top=coords_1_mid_top(:,1); 
Yi_1_mid_top=coords_1_mid_top(:,2);
Vi_1_mid_top = gIDW(Xc_1_mid_top,Yc_1_mid_top,Vc_1_mid_top,Xi_1_mid_top,Yi_1_mid_top,-2);

%Using the coordinates from the line drawn above, input the interpolated
%topography into the large-scale topography variable.
for i=1:size(coords_1_mid_top,1)
    for j=coords_1_mid_top(i,1)
        for k=coords_1_mid_top(i,2)
            topo(j,k)=Vi_1_mid_top(i);
            break
        end
    end
end


% elevation between mid point and input_topo_2, top
% Create variables with rows/cols of top of input_topo_2 (right edge) and mid point
% and elevation at these locations
Xc_1_2_top = [mid_point_1_top_row,rows_2_top_TOPO(end)];
Yc_1_2_top = [mid_point_1_top_col,cols_2_top_TOPO(end)];
Vc_1_2_top = [mid_point_1_top_elev, avg_elev_2_top];

% Create a line with 800 points using the coordinates from above.
[fill_rows_1_2_top, fill_cols_1_2_top, I_1_2_top]=improfile(topo,Xc_1_2_top,Yc_1_2_top,800);

% Use the coordinates from this line along with the coordinates and
% elevations of the top endpoints of input_topo_2 and midpoint to calculate
% elevations along this line. We use the gIDW which computes the inverse
% weighted interpolation. 
coords_1_2_top=round([fill_rows_1_2_top, fill_cols_1_2_top]);
Xi_1_2_top=coords_1_2_top(:,1); 
Yi_1_2_top=coords_1_2_top(:,2);
Vi_1_2_top = gIDW(Xc_1_2_top,Yc_1_2_top,Vc_1_2_top,Xi_1_2_top,Yi_1_2_top,-2);

%Using the coordinates from the line drawn above, input the interpolated
%topography into the large-scale topography variable.
for i=1:size(coords_1_2_top,1)
    for j=coords_1_2_top(i,1)
        for k=coords_1_2_top(i,2)
            topo(j,k)=Vi_1_2_top(i); 
            break
        end
    end
end

%%% SAME STEPS FROM TOPO 1 TO MIDPOINT TOP %%% 
% repeating the steps from above, now for the bottom edge between the
% input_topo_1 and input_topo_2
% elevation half way between input_topo_1 and 2, bottom
mid_point_1_bot_row = round((rows_1_bot_TOPO(1)+rows_2_bot_TOPO(end))/2);
mid_point_1_bot_col = round((cols_1_bot_TOPO(1)+cols_2_bot_TOPO(end))/2);
mid_point_1_bot_elev = (avg_elev_1_bot+avg_elev_2_bot)/2;
topo(mid_point_1_bot_row,mid_point_1_bot_col)=mid_point_1_bot_elev;

% elevation between input_topo_1 and midpoint, bottom
Xc_1_mid_bot = [rows_1_bot_TOPO(1),mid_point_1_bot_row];
Yc_1_mid_bot = [cols_1_bot_TOPO(1),mid_point_1_bot_col ];
Vc_1_mid_bot = [avg_elev_1_bot, mid_point_1_bot_elev];
[fill_rows_1_mid_bot, fill_cols_1_mid_bot, I_1_mid_bot]=improfile(topo,Xc_1_mid_bot,Yc_1_mid_bot,800);
coords_1_mid_bot=round([fill_rows_1_mid_bot, fill_cols_1_mid_bot]);
Xi_1_mid_bot=coords_1_mid_bot(:,1); 
Yi_1_mid_bot=coords_1_mid_bot(:,2);
Vi_1_mid_bot = gIDW(Xc_1_mid_bot,Yc_1_mid_bot,Vc_1_mid_bot,Xi_1_mid_bot,Yi_1_mid_bot,-2);

for i=1:size(coords_1_mid_bot,1)
    for j=coords_1_mid_bot(i,1)
        for k=coords_1_mid_bot(i,2)
            topo(j,k)=Vi_1_mid_bot(i); 
            break
        end
    end
end

% elevation beteween mid point and input_topo_2, bottom
Xc_1_2_bot = [mid_point_1_bot_row,rows_2_bot_TOPO(end)];
Yc_1_2_bot = [mid_point_1_bot_col,cols_2_bot_TOPO(end)];
Vc_1_2_bot = [mid_point_1_bot_elev, avg_elev_2_bot];
[fill_rows_1_2_bot, fill_cols_1_2_bot, I_1_2_bot]=improfile(topo,Xc_1_2_bot,Yc_1_2_bot,800);
coords_1_2_bot=round([fill_rows_1_2_bot, fill_cols_1_2_bot]);
Xi_1_2_bot=coords_1_2_bot(:,1); 
Yi_1_2_bot=coords_1_2_bot(:,2);
Vi_1_2_bot = gIDW(Xc_1_2_bot,Yc_1_2_bot,Vc_1_2_bot,Xi_1_2_bot,Yi_1_2_bot,-2);

for i=1:size(coords_1_2_bot,1)
    for j=coords_1_2_bot(i,1)
        for k=coords_1_2_bot(i,2)
            topo(j,k)=Vi_1_2_bot(i); 
            break
        end
    end
end

%% Adding average elevations to points between Landlab topographies, ridgelines

%%% SAME STEPS FROM TOPO 1 TO MIDPOINT TOP %%% 
% repeat the steps used to draw the top and bottom edges between
% input_topo_1 and input_topo_2 but now for the ridgeline between them.

% NOTE: for the ridgeline we add a 1/4 point and 3/4 point along with the
% midpoint. This allows for a smoother gradient in elevations between the
% two domains. The steps for adding lines between the points remains the
% same as those done above. 

% Elevation between input_topo_1 and input_topo_2, ridgeline
mid_point_1_row=round((rows_1_L_TOPO(round(length(rows_1_L_TOPO)/2))+rows_2_R_TOPO(round(length(rows_2_R_TOPO)/2)))/2);
mid_point_1_col=round((cols_1_L_TOPO(round(length(cols_1_L_TOPO)/2))+cols_2_R_TOPO(round(length(cols_2_R_TOPO)/2)))/2);
mid_point_1_elev=(avg_elev_1_L+avg_elev_2_R)/2;
topo(mid_point_1_row, mid_point_1_col)=mid_point_1_elev;
% find the coordinate for the mid point by finding the middle row/col of the
% left side of input_topo_1 and adding it to the middle row/col of the right side of
% input_topo_2 and dividing this by 2

qu_point_1_row=round((rows_1_L_TOPO(round(length(rows_1_L_TOPO)/2))+mid_point_1_row)/2);
qu_point_1_col=round((cols_1_L_TOPO(round(length(cols_1_L_TOPO)/2))+mid_point_1_col)/2);
qu_point_1_elev=(avg_elev_1_L+(avg_elev_1_L+avg_elev_2_R)/2)/2;
topo(qu_point_1_row,qu_point_1_col)=qu_point_1_elev;
% find the coordinate of the quarter point by finding the middle row/col of the
% left side of input_topo_1 adding it to the mid point point coordinate and
% dividing by 2

tqu_point_1_row=round((rows_2_R_TOPO(round(length(rows_2_R_TOPO)/2))+mid_point_1_row)/2);
tqu_point_1_col=round((cols_2_R_TOPO(round(length(cols_2_R_TOPO)/2))+mid_point_1_col)/2);
tqu_point_1_elev=(avg_elev_2_L+(avg_elev_1_L+avg_elev_2_R)/2)/2;
topo(tqu_point_1_row, tqu_point_1_col)=tqu_point_1_elev;

% Filling in ridegline from input_topo_1 to 1/4 point    
Xc_1_qu=double([rows_1_L_TOPO(round(length(rows_1_L_TOPO)/2)),qu_point_1_row]); % location of rows of left point and right point
Yc_1_qu=double([cols_1_L_TOPO(round(length(cols_1_L_TOPO)/2)), qu_point_1_col]); % location of cols of left point and right point
Vc_1_qu=[avg_elev_1_L, qu_point_1_elev]; % topography values at left point and right point
[fill_rows_1_qu, fill_cols_1_qu, I_1_qu]=improfile(topo,Xc_1_qu,Yc_1_qu,800); % locations of cells to be filled in
coords_1_qu=round([fill_rows_1_qu, fill_cols_1_qu]); % define matrix of coordinates to be filled
Xi_1_qu=(double(coords_1_qu(:,1))); % rows to fill
Yi_1_qu=(double(coords_1_qu(:,2))); % cols to fill
Vi_1_qu = gIDW(Xc_1_qu,Yc_1_qu,Vc_1_qu,Xi_1_qu,Yi_1_qu,-2); % inverse distance weighting between left and right points 

for i=1:size(coords_1_qu,1)
    for j=coords_1_qu(i,1)
        for k=coords_1_qu(i,2)
            topo(j,k)=Vi_1_qu(i); %filling in topo with new values
            break
        end
    end
end


% Filling in ridgeline from 1/4 point to 1/2 point
Xc_1_qu_mid=double([qu_point_1_row, mid_point_1_row]);
Yc_1_qu_mid=double([qu_point_1_col, mid_point_1_col]);
Vc_1_qu_mid=[qu_point_1_elev, mid_point_1_elev];
[fill_rows_1_qu_mid, fill_1_cols_qu_mid, I_1_qu_mid]=improfile(topo, Xc_1_qu_mid,Yc_1_qu_mid,800);
coords_1_qu_mid=round([fill_rows_1_qu_mid,fill_1_cols_qu_mid]);
Xi_1_qu_mid=double(coords_1_qu_mid(:,1));
Yi_1_qu_mid=double(coords_1_qu_mid(:,2));
Vi_1_qu_mid = gIDW(Xc_1_qu_mid, Yc_1_qu_mid,Vc_1_qu_mid,Xi_1_qu_mid,Yi_1_qu_mid,-2);

for i=1:size(coords_1_qu_mid,1)
    for j=coords_1_qu_mid(i,1)
        for k=coords_1_qu_mid(i,2)
            topo(j,k)=Vi_1_qu_mid(i);
            break
        end
    end
end

% Filling in ridgeline from 1/2 point to 3/4 point
Xc_1_mid_tqu=double([mid_point_1_row,tqu_point_1_row]);
Yc_1_mid_tqu=double([mid_point_1_col,tqu_point_1_col]);
Vc_1_mid_tqu=[mid_point_1_elev, tqu_point_1_elev];
[fill_rows_1_mid_tqu, fill_cols_1_mid_tqu, I_1_mid_tqu]=improfile(topo, Xc_1_mid_tqu, Yc_1_mid_tqu,800);
coords_1_mid_tqu=round([fill_rows_1_mid_tqu, fill_cols_1_mid_tqu]);
Xi_1_mid_tqu=double(coords_1_mid_tqu(:,1));
Yi_1_mid_tqu=double(coords_1_qu_mid(:,2));
Vi_1_mid_tqu = gIDW(Xc_1_mid_tqu, Yc_1_mid_tqu,Vc_1_mid_tqu,Xi_1_mid_tqu,Yi_1_mid_tqu,-2);

for i=1:size(coords_1_mid_tqu,1)
    for j=coords_1_mid_tqu(i,1)
        for k=coords_1_mid_tqu(i,2)
            topo(j,k)=Vi_1_mid_tqu(i);
            break
        end
    end
end

% Filling in rideline from 3/4 point to input_topo_2
Xc_1_tqu_2=double([tqu_point_1_row,rows_2_R_TOPO(round(length(rows_2_R_TOPO)/2))]);
Yc_1_tqu_2=double([tqu_point_1_col,cols_2_R_TOPO(round(length(cols_2_R_TOPO)/2))]);
Vc_1_tqu_2=[tqu_point_1_elev, avg_elev_2_R];
[fill_rows_1_tqu_2, fill_cols_1_tqu_2, I_1_tqu_2]=improfile(topo, Xc_1_tqu_2, Yc_1_tqu_2,800);
coords_1_tqu_2=round([fill_rows_1_tqu_2, fill_cols_1_tqu_2]);
Xi_1_tqu_2=double(coords_1_tqu_2(:,1));
Yi_1_tqu_2=double(coords_1_tqu_2(:,2));
Vi_1_tqu_2 = gIDW(Xc_1_tqu_2, Yc_1_tqu_2,Vc_1_tqu_2,Xi_1_tqu_2,Yi_1_tqu_2,-2);

for i=1:size(coords_1_tqu_2,1)
    for j=coords_1_tqu_2(i,1)
        for k=coords_1_tqu_2(i,2)
            topo(j,k)=Vi_1_tqu_2(i);
            break
        end
    end
end

%% Interpolating topography between Landlab domains
% Here we use the lines of coordinates and elevations created above and interpolate
% them on a mesh that is then added to the large-scale topography file. 

% Interpolate topography between input_topo_1 and input_topo_2
% Create a matrix of all the known coordinates (rows, columns) and elevations

% concacate the rows for the top/ridgeline/bottom  fill topographies
XZ_1 = [coords_1_mid_top(:,1);coords_1_2_top(:,1);coords_1_qu(:,1);coords_1_qu_mid(:,1);coords_1_mid_tqu(:,1);coords_1_tqu_2(:,1);coords_1_mid_bot(:,1);coords_1_2_bot(:,1)]; 

% concacate the cols for the top/ridgeline/bottom fill topographies
YZ_1 = [coords_1_mid_top(:,2);coords_1_2_top(:,2);coords_1_qu(:,2);coords_1_qu_mid(:,2);coords_1_mid_tqu(:,2);coords_1_tqu_2(:,2);coords_1_mid_bot(:,2);coords_1_2_bot(:,2)]; 

% concacate fill topography values
ZZ_1 = [Vi_1_mid_top;Vi_1_2_top;Vi_1_qu;Vi_1_qu_mid;Vi_1_mid_tqu;Vi_1_tqu_2;Vi_1_mid_bot;Vi_1_2_bot]; 

% Create the interpolant from the known data
Fz_1 = scatteredInterpolant(XZ_1,YZ_1,ZZ_1); 

% Build a mesh of query points that represent the space between input_topo_1 and input_topo_2

%repeat row coordinates from large-scale topography across the number of columns between input_topo_1 and north end point
XZ_mesh_1 = (repmat(((rows_1_L_TOPO(1):rows_2_R_TOPO(end))),3554,1))';

%repeat column coordiantes from large-scale topography across the number of rows between input_topo_1 and north end point
YZ_mesh_1 = (repmat(((cols_2_R_TOPO(1):cols_1_L_TOPO(end))),1408,1));

%Evaulate the interpolant at these query points 
Vq_1 = Fz_1(XZ_mesh_1, YZ_mesh_1); 

% add the interpolated mesh to the large-scale topography variable
for i = 1:size(Vq_1,1)
    for j = 1:size(Vq_1,2)
        if Vq_1(i,j)<1^-3 % at the edges of the interpolated mesh, set very small values to 0
            Vq_1(i,j)=0;
        end
        % if the topography variable is 0 add the interpolated topography
        if topo(i+5506,j+26696)==0 % NOTE: the location that the interpolated topography is inserted must be manually selected
            topo(i+5506,j+26696)=Vq_1(i,j);
        % if the topography variable is not 0, you have reached input_topo_1/2 and
        % should not overwrite the Landlab topography 
        elseif topo(i+5506,j+26696)>0
            
        end
    end
end

%% Write out the final topography file as a netcdf that will be used in WRF %% 
rows=size(topo,1);
cols=size(topo,2);
nccreate('geo_em.coor_interpolated_topo.nc','HGT_M','Dimensions',{'rows',rows,'cols',cols})
ncwrite('geo_em.coor_interpolated_topo.nc','HGT_M', topo)