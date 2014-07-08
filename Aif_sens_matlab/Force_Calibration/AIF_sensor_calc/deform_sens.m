
function [dx_ends, def_centers,Cal_im_back,Def_im_back] = deform_sens(cal_image_matrix, cal_image_info, deform_image_info, tag_edge_length, Tags_of_interest)

%SHOULD DO A DIMENSION CHECK TO ENSURE ALL TAGS ARE VISIBLE

%I expect Cal_image and deform_image to have the form:
%  [tagA tagB tagC;...
%   pcA  pcB  pcB ;...
%   p1A  p1B  p1C ;...
%    ............
%   p4A  p4B  p4C ];

%Calibration image needs to be taken with a rule in the background,
%calibration image needs rule of 11in across the top of the frame, so 
%the sensor can be placed on a sheet of paper (landscape) and the 11in side
%of the paper should go across the top of the camera. And i assume the
%camera is level.

% Tag information
% Tags_of_interest = [2;...   %tag center
%                     98;...  %tag1
%                     62;...  %tag2
%                     50;...  %tag3
%                     38];    %tag4


%Let Cal_dist be arranged as:
% [dist_cent_tag1_mm, dist_cent_tag2_mm, dist_cent_tag3_mm,... 
%  dist_cent_tag4_mm, dist_cent_tag5_mm, dist_cent_tag6_mm];
%


%tag_edge_length: length of tag edge (black edge) in mm 

%% Analyize calibration image first



  
Cal_image_arranged = [];
tags_ordered = [];
for i = 1:length(Tags_of_interest)   
[~, tag] = find(cal_image_info(1,:) == Tags_of_interest(i)); %this assumes tags of 
tags_ordered(end+1) = tag;
end

for i = 1:length(Tags_of_interest)
Cal_image_arranged(:,end+1) = cal_image_info(:,tags_ordered(i));
end


Cal_im_back = Cal_image_arranged;

Cal_image_arranged = double(Cal_image_arranged);

 
        
        

%% Rearrange deformed image tags
   

Def_image_arranged = [];
tags_ordered = [];
for i = 1:length(Tags_of_interest)   
[~, tag] = find(deform_image_info(1,:) == Tags_of_interest(i)); %this assumes tags of 
tags_ordered(end+1) = tag;
end

for i = 1:length(Tags_of_interest)
Def_image_arranged(:,end+1) = deform_image_info(:,tags_ordered(i));
end
Def_im_back = Def_image_arranged;
                  
Def_image_arranged = double(Def_image_arranged);






%% Now find the skew matrix between the center tags

%Eqn: Deflect_cent = transform* calibrate_cent
%     transform = pinv(calibrate_cent)*Deflect_cent


%as it is now these are vectors from the corner of the screen, i need from
%a set point in the sensor, which i choose to be the p1 of the center

Cal_sensor_frame_centertag = Cal_image_arranged(4:end,1)-repmat(Cal_image_arranged(4:5,1),size(Cal_image_arranged(4:end,1),1)/2,size(Cal_image_arranged(4:end,1),2));
Def_sensor_frame_centertag = Def_image_arranged(4:end,1)-repmat(Def_image_arranged(4:5,1),size(Def_image_arranged(4:end,1),1)/2,size(Def_image_arranged(4:end,1),2));

Cal_sensor_frame_centertag_trans = reshape(Cal_sensor_frame_centertag,2,4);

Def_sensor_frame_centertag_trans = reshape(Def_sensor_frame_centertag,2,4);

S_mat = Def_sensor_frame_centertag_trans*pinv(Cal_sensor_frame_centertag_trans); %transform of center tag





%% Find projected position of all tag corners (p1 (not pc))



Cal_sensor_frame = Cal_image_arranged(4:end,1:end)-repmat(Cal_image_arranged(4:5,1),size(Cal_image_arranged(4:end,1:end),1)/2,size(Cal_image_arranged(4:end,1:end),2));
Def_sensor_frame = Def_image_arranged(4:end,1:end)-repmat(Def_image_arranged(4:5,1),size(Def_image_arranged(4:end,1:end),1)/2,size(Def_image_arranged(4:end,1:end),2));
 

Cal_sensor_tag_corners = Cal_sensor_frame(1:2,1:end);

Def_sensor_tag_corners = Def_sensor_frame(1:2,1:end);

Projected_sensor_tag_corners = S_mat * Cal_sensor_tag_corners;

%this is the linear transformation: [dx; dy] of the reference frames (each tag)
T = Def_sensor_tag_corners - Projected_sensor_tag_corners; 




%now convert all other corners (other than primary corners) to projected
%plane for comparison

Cal_end_ref_frames = Cal_sensor_frame(3:end,:) - repmat(Cal_sensor_tag_corners,3,1);



Def_end_ref_frames = Def_sensor_frame(3:end,:) - repmat(Def_sensor_tag_corners,3,1);


%Find rotation matrix between cal and deformed 
% Every cell of R will contain the rotation matrix for that end effector
R = cellmat;

col_cerf = (size(Cal_end_ref_frames,1)*size(Cal_end_ref_frames,2))/2;
Cal_end_ref_frames = reshape(Cal_end_ref_frames,2,col_cerf);

col_derf = (size(Def_end_ref_frames,1)*size(Def_end_ref_frames,2))/2;
Def_end_ref_frames = reshape(Def_end_ref_frames,2,col_derf);

for i = 1: (size(Cal_end_ref_frames,2))/3
   R{i} =   (S_mat*Cal_end_ref_frames(:,(1:3)*i)) * pinv(Def_end_ref_frames(:,(1:3)*i));
end

%So now at this point I have R and T


%% Find Difference

a_cent_def_possible = Def_image_arranged(2:3,:) - Def_image_arranged(4:5,:);


deflection = [];

I = [1 0;...
     0 1];

for i = 1:length(R)

deflection(:,end+1) = [ (R{i} - I),T(:,i)]*[a_cent_def_possible(:,i); 1];
end


%% Transform deflection back into calibrated frame 

def_centers = pinv(S_mat)*deflection;

%% Find true deflection in mm

scale = mean([norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(6:7,1)]),...
             norm([Cal_image_arranged(6:7,1)-Cal_image_arranged(8:9,1)]),...
             norm([Cal_image_arranged(8:9,1)-Cal_image_arranged(10:11,1)]),...
             norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(10:11,1)])]);



%Convert this distance to metric units
%image matrix is  [y x] = size(I), and starts from upper left and corner

%  dx_ends = int64((tag_edge_length/scale).*def_centers);  %this is distance in mm (if scale is 11in)
  
dx_ends = (tag_edge_length/scale).*def_centers;
 



end