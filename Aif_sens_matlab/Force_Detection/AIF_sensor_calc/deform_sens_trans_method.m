%Deformation Calculation using Peter Corke (camera transformation) Method

function [dx_ends, def_centers,Cal_im_back,Def_im_back] = deform_sens_trans_method(cal_image_matrix, cal_image_info, deform_image_info, tag_edge_length, Tags_of_interest, Calib_matrix)



%I expect Cal_image and deform_image to have the form:
%  [tagA tagB tagC;...
%   pcA  pcB  pcB ;...  %position of center of each tag
%   p1A  p1B  p1C ;...  %position of corner 1 of each tag
%    ............
%   p4A  p4B  p4C ];

% Tag information

%tag 1: 13868794921
%tag 2: 38943287133
%tag 3: 32013986588
%tag 4: 37856557532
%tag 5: 47628935088
%tag 6: 7661251159
%tag center: 59366215950

%Let Cal_dist be arranged as:
% [dist_cent_tag1_mm, dist_cent_tag2_mm, dist_cent_tag3_mm,...
%  dist_cent_tag4_mm, dist_cent_tag5_mm, dist_cent_tag6_mm];
%


%tag_edge_length: length of tag edge (black edge) in mm

%% Analyize calibration image first




Cal_image_arranged = [];
tags_ordered = [];
for i = 1:length(Tags_of_interest)
    [~, tag] = find(cal_image_info == Tags_of_interest(i)); %this assumes tags of
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
    [~, tag] = find(deform_image_info == Tags_of_interest(i)); %this assumes tags of
    tags_ordered(end+1) = tag;
end

for i = 1:length(Tags_of_interest)
    Def_image_arranged(:,end+1) = deform_image_info(:,tags_ordered(i));
end
Def_im_back = Def_image_arranged;

Def_image_arranged = double(Def_image_arranged);







%% Camera Transformation method
%In this method, the center tag in the calibration image is used to find
%"pose" of the camera, this "pose" is then used to establish calibrated
%position of end effectors. Then new "pose" of camera is found in deformed
%image, and based on saved position of end effectors from calibration, the
%difference of end effectors is found.

%ground assignment for center tag: (TL := taglength)
% _____| X  | Y
% Pc   |TL/2|TL/2
% P1   | 0  | 0
% P2   |TL  | 0
% P3   |TL  |TL
% P4   | 0  |TL
%
% pc_cent_cal = Cal_image_arranged(2:3,1); %2x1 [x,y] of center
% p1_cent_cal = Cal_image_arranged(4:5,1);
% p2_cent_cal = Cal_image_arranged(6:7,1);
% p3_cent_cal = Cal_image_arranged(8:9,1);
% p4_cent_cal = Cal_image_arranged(10:11,1);



TL = tag_edge_length;

%cent_location = [TL/2, TL/2;...
%                  0   ,    0;...
%                  TL  ,    0;...
%                  TL  ,   TL;...
%                  0   ,   TL];  %locations of center tag corners and center
%
cent_location = [0   ,    0;...
    -TL/2,-TL/2;...
    TL/2,-TL/2;...
    TL/2, TL/2;...
    -TL/2, TL/2];  %locations of center tag corners and center


%Undistort Cal_image_arranged
Cal_image_stretch = reshape(Cal_image_arranged(2:end,:),2,(numel(Cal_image_arranged)-size(Cal_image_arranged,2))/2);

Cal_image_stretch = [Cal_image_stretch;...
    ones(1,size(Cal_image_stretch,2))];

Cal_img_undistort = inv(Calib_matrix)*Cal_image_stretch;

Cal_img_undistort = Cal_img_undistort(1:2,:);



M = [];
for i = 1:5
    
    x = Cal_img_undistort(1,i);
    y = Cal_img_undistort(2,i);
    X = cent_location(i,1);
    Y = cent_location(i,2);
    
    M = [M;...
        -X, 0, x*X, -Y, 0, x*Y, -1, 0, x;...
        0, -X, y*X, 0, -Y, y*Y, 0, -1, y];
end

[S, U, V] = svd(M);

A = V(:,end); %this is the eigenvector of the nullspace that best approximates rotation and translation

if A(end)<0
    A = -A;
end

%normalize
% r_norm = norm(A(1:3,1));
% A = (1/r_norm)*A;

%transform (this assumes all tags are in the plane)
A_cal = reshape(A,3,3); %this transform basically prescribes the position of the camera with respect to the center tag


%Now find Position of all

% Cal_image_stretch = reshape(Cal_image_arranged(2:end,:),2,(numel(Cal_image_arranged)-size(Cal_image_arranged,2))/2);

% Cal_image_stretch = [Cal_image_stretch;...
%     ones(1,size(Cal_image_stretch,2))];

Cal_image_world_frame_stretch = pinv(A_cal)*[Cal_img_undistort;ones(1,size(Cal_img_undistort,2))]; %Now we have all the points of every tag in the world frame (since we calibrated with the center tag)

Cal_image_world_frame_stretch = Cal_image_world_frame_stretch(1:2,:); %delete the row of ones that are non-essential

Cal_image_world_frame = reshape(Cal_image_world_frame_stretch, size(Cal_image_arranged(2:end,:),1),size(Cal_image_arranged(2:end,:),2)); %now put back into an array: only this array is in world frame for every point


%% Now look at the deformed image


%Undistort Def_image_arranged
Def_image_stretch = reshape(Def_image_arranged(2:end,:),2,(numel(Def_image_arranged)-size(Def_image_arranged,2))/2);

Def_image_stretch = [Def_image_stretch;...
    ones(1,size(Def_image_stretch,2))];

Def_img_undistort = inv(Calib_matrix)*Def_image_stretch;

Def_img_undistort = Def_img_undistort(1:2,:);


M = [];

for i = 1:5
    
    x = Def_img_undistort(1,i);
    y = Def_img_undistort(2,i);
    X = cent_location(i,1);
    Y = cent_location(i,2);
    
    M = [M;...
        -X, 0, x*X, -Y, 0, x*Y, -1, 0, x;...
        0, -X, y*X, 0, -Y, y*Y, 0, -1, y];
end

[S, U, V] = svd(M);

A = V(:,end); %this is the eigenvector of the nullspace that best approximates rotation and translation

if A(end)<0
    A = -A;  %this step is done becuase there is an equivalent solution behind the camera
end

%normalize
% r_norm = norm(A(1:3,1));
% A = (1/r_norm)*A;

%transform (this assumes all tags are in the plane)
A_def = reshape(A,3,3); %this transform basically prescribes the position of the camera with respect to the center tag

%Now find Position of all



Def_image_world_frame_stretch = pinv(A_def)*[Def_img_undistort;ones(1,size(Def_img_undistort,2))]; %Now we have all the points of every tag in the world frame (since we calibrated with the center tag)

Def_image_world_frame_stretch = Def_image_world_frame_stretch(1:2,:); %delete the row of ones that are non-essential

Def_image_world_frame = reshape(Def_image_world_frame_stretch, size(Def_image_arranged(2:end,:),1),size(Def_image_arranged(2:end,:),2)); %now put back into an array: only this array is in world frame for every point





%% Now find the deflection and rotation


Deflection_array = Def_image_world_frame - Cal_image_world_frame;

%Note: the deflections of every tag can be found from top 2 rows, the
%rotation can be determined from rest


dx_ends = Deflection_array(1:2,:);
% 
% diff_uv_cal = Calib_matrix*A_cal*[dx_ends; ones(1,size(dx_ends,2))];
% 
% 
% def_check = Calib_matrix*A_cal*[Cal_image_world_frame(1:2,:);ones(1,size(Cal_image_world_frame,2))];

%now find the pixel deflection for display:

scale = mean([norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(6:7,1)]),...
             norm([Cal_image_arranged(6:7,1)-Cal_image_arranged(8:9,1)]),...
             norm([Cal_image_arranged(8:9,1)-Cal_image_arranged(10:11,1)]),...
             norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(10:11,1)])]);

         def_centers_old = (scale/TL).*[dx_ends(1,:);dx_ends(2,:)]; %needs flipping for display

         
         
 %% need to use the method below, and adjust in the print       
         

cal_img_def_pixels = A_cal* [Def_image_world_frame(1:2,:);ones(1,size(Def_image_world_frame,2))]; %deformed position from calibration perspective (to portray on cal image)


def_centers = Calib_matrix*cal_img_def_pixels;

def_in_calib_frame =  Calib_matrix * A_cal * inv(A_def) * inv(Calib_matrix) * Def_image_stretch;

def_in_calib_frame = def_in_calib_frame(1:2,:);

def_centers = reshape(def_in_calib_frame, size(Def_image_arranged(2:end,:),1),size(Def_image_arranged(2:end,:),2));
def_centers = def_centers(1:2,:);

%   def_centers = def_centers(1:2,:)-Cal_image_arranged(2:3,:);

%  def_centers = Calib_matrix*A_cal*[dx_ends;]













