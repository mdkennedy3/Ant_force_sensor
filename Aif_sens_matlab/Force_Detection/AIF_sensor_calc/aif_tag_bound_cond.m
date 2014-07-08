function [tag_ang_BC, Ls_mat, tag_init_position, tag_pos_bound_cond] = aif_tag_bound_cond(cal_image_matrix, cal_image_info, deform_image_info, tag_edge_length, Tags_of_interest, dx_ends)



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

%% Step 1: Analyize & Rearrange calibration image, and deformed image tags

%%%%%%%%%%Cal_image%%%%%%%%%%%%%%%%%%%%  
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

%%%%%%%%%%Def_image%%%%%%%%%%%%%%%%%%%%
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



%% Step 2: Turn tags into vectors, where rotation is assumed to be about center of tag

%%%%Calibration Image%%%%%%%%%%%%%%%%%%%%
%Note: Each tag will have four vectors, (from the center to each corner),
%where vectors point from the center to the corners

% Cal_tag_vectors =
% center     T1   .... T4
% [ (p1-pc) ....
% [ (p2-pc) ....
% [ (p3-pc) ....
% [ (p4-pc) ....

Cal_img_vecs = [];

for i = 1:size(Cal_image_arranged,2)
    
cent = Cal_image_arranged(2:3,i);
p1 = Cal_image_arranged(4:5,i);
p2 = Cal_image_arranged(6:7,i);
p3 = Cal_image_arranged(8:9,i);
p4 = Cal_image_arranged(10:11,i);

    
Cal_img_vecs(:,i) = [(p1-cent)/norm((p1-cent),2);...
                     (p2-cent)/norm((p2-cent),2);...
                     (p3-cent)/norm((p3-cent),2);...
                     (p4-cent)/norm((p4-cent),2)];
end


%%%%%%%%%%%Deformed Image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Def_img_vecs = [];

for i = 1:size(Def_image_arranged,2)
    
cent = Def_image_arranged(2:3,i);
p1 = Def_image_arranged(4:5,i);
p2 = Def_image_arranged(6:7,i);
p3 = Def_image_arranged(8:9,i);
p4 = Def_image_arranged(10:11,i);

    
Def_img_vecs(:,i) = [(p1-cent)/norm((p1-cent),2);...
                     (p2-cent)/norm((p2-cent),2);...
                     (p3-cent)/norm((p3-cent),2);...
                     (p4-cent)/norm((p4-cent),2)];
end



%% Step 3: Find average rotation angles for each tag  
% Below I use quaternions to find the angle of rotation of each tag
% assuming that they are each rotated about their center through a vector
% that is on the z axis (vector normal to tag), quaternion: q = q0 + q3(k)
% hence q = cos(th/2) + sin(th/2) (k)

Tag_thetas = [];

for i = 1:size(Cal_image_arranged,2)
   
    th_vec = [];
    for j = 1:4
    v1 = Cal_img_vecs(j*2-1,i);
    v2 = Cal_img_vecs(j*2,i);
    w1 = Def_img_vecs(j*2-1,i);
    w2 = Def_img_vecs(j*2,i);
   
    th_vec(end+1) = 2*acos(sqrt( (w2 + w1*(v1/v2))/(2*((v1^2)/(v2)+ v2)) + 0.5 )) ; %rotation from Calibrated Frame to Deformed frame
    end
    
    Tag_thetas(i,1) = mean(th_vec); 
    
end


%% Step 4: Find respective angles for each tag T1, ..., T4 with respect to rotation of Tc (center tag)


tag_ang_BC = zeros(size(Cal_image_arranged,2),1);

%center tag rotation is stationary as it is used for reference (rotating
%hanging weight back into cal_image frame)
tag_ang_BC(1) = Tag_thetas(1);

tag_ang_BC(2:end) = Tag_thetas(2:end) - repmat(Tag_thetas(1),length(Tag_thetas)-1,1); 




%% Get L* (at rest length of every coil)


scale = mean([norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(6:7,1)]),...
             norm([Cal_image_arranged(6:7,1)-Cal_image_arranged(8:9,1)]),...
             norm([Cal_image_arranged(8:9,1)-Cal_image_arranged(10:11,1)]),...
             norm([Cal_image_arranged(4:5,1)-Cal_image_arranged(10:11,1)])]);

Ls_mat = [];         

for i = 2:size(Cal_image_arranged,2)
Ls_mat(end+1,1) = (tag_edge_length/scale).*abs([norm((Cal_image_arranged(2:3,i)-Cal_image_arranged(2:3,1)),2)]); %approximation of length of coil in mm
end



%% Get Position Boundary Condidtions
%it is assumed that the B.C. for the center tag is: dy/dx(0) = 0,, y(0) = 0


%find L* with direction (vector), then add deflection to that

tag_init_position = []; %initial position of tags
for i = 2:size(Cal_image_arranged,2)
    tag_init_position(1:2,end+1) = (tag_edge_length/scale).*[(Cal_image_arranged(2:3,i)-Cal_image_arranged(2:3,1))]; %approximation of length of coil in mm
end

tag_pos_bound_cond = [];

for i = 2:size(Cal_image_arranged,2)
   tag_pos_bound_cond(1:2,end+1) = tag_init_position(1:2,i-1) + dx_ends(1:2,i); %Boundary conditions [x(L)_T1, x(L)_T2...
                                                                           %                     y(L)_T1, y(L)_T2...] 
end

%Now all of these are distances from center tag, and don't worry about
%orientation, because forces due to gravity will be rotated into the
%calibrated frame, so when a tag of interest is considered then the
%rotation will be consistent with the physics in the deformed frame.



