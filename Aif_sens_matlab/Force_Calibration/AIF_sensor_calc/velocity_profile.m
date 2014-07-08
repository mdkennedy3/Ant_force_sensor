function [vel_prof] = velocity_profile(def_img_info_init_frame,F_init, def_img_info_final_frame,F_final, tag_edge_length, Tags_of_interest_init,Tags_of_interest_final, FPS)

%Input:
%
% img_info
%  [tagA tagB tagC;...
%   pcA  pcB  pcB ;...
%   p1A  p1B  p1C ;...
%    ............
%   p4A  p4B  p4C ];
%
%
%
% def_img_info_init_frame: Image info from deformed 
% F_init: frame number of 1st frame
% 
% def_img_info_final_frame:
% F_final: frame number of second frame
%
% tag_edge_length: length in mm of edge of tag
% FPS: frame per second (29.97 most likely)

%Output: 

%
% vel_prof = [tagA      tagB        tagC;...  %tag names
%             vel_c_A   vel_c_B     vel_c_C;... %velocities of centers
%             w_vel_A   w_vel_B     w_vel_C;... %angular velocities of each
%             frame#     frame#       frame#] %NOTE these are the same
%             (note: given frame N and N+1 that frame# is FrameN (as its
%             vel vector from the inital point)
%

%% Order the tags appropriately 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%handle initial frame%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Def_image_init_arranged = [];
tags_ordered = [];
for i = 1:length(Tags_of_interest_init)   
[~, tag] = find(def_img_info_init_frame == Tags_of_interest_init(i)); %this assumes tags of 
tags_ordered(end+1) = tag;
end

for i = 1:length(Tags_of_interest_init)
Def_image_init_arranged(:,end+1) = def_img_info_init_frame(:,tags_ordered(i));
end
                  
Def_image_init_arranged = double(Def_image_init_arranged);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Rearrange deformed image tags
%%%%%%%%%%%%%%%%%%%%%handle final frame%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Def_image_final_arranged = [];
tags_ordered = [];
for i = 1:length(Tags_of_interest_final)   
[~, tag] = find(def_img_info_final_frame == Tags_of_interest_final(i)); %this assumes tags of 
tags_ordered(end+1) = tag;
end

for i = 1:length(Tags_of_interest_final)
Def_image_final_arranged(:,end+1) = def_img_info_final_frame(:,tags_ordered(i));
end
                  
Def_image_final_arranged = double(Def_image_final_arranged);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Now truncate so that they only compare the same tags

if (length(Tags_of_interest_final) < length(Tags_of_interest_init))
    
    Def_image_init_arranged = Def_image_init_arranged(:,1:length(Tags_of_interest_final));
    
elseif (length(Tags_of_interest_final) > length(Tags_of_interest_init))
        
    Def_image_final_arranged = Def_image_final_arranged(:,1:length(Tags_of_interest_init)); 
    
end
    
   
%% 1. Find the conversion btw pixel length and tag length (mm)


%realize that what I've done is found the distance from every corner to
%center, then realizing that this average distance is equal to
%(1/sqrt(2))*edge_length as tags are always square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corn_to_center_pix = Def_image_final_arranged(4:end,:) - repmat(Def_image_final_arranged(2:3,:),size(Def_image_final_arranged(4:end,:),1)/2,1);
%now find distances to each (as this is x,y distance) -> dist = sqrt(x^2 +
%y^2)
corn_to_center_pix_stretch= reshape(corn_to_center_pix,2, numel(corn_to_center_pix)/2);
%find l2 distance
corn_to_center_pix_dist = sqrt(sum((corn_to_center_pix_stretch.*corn_to_center_pix_stretch),1));

%then reshape and find the mean distance from corner to center for every
%tag
corn_to_center_pix_dist = mean(reshape(corn_to_center_pix_dist,4,int8(numel(corn_to_center_pix_dist)/4)),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dimension of below is 1xn where n is number of tags
pixel_mm_conv = (tag_edge_length*(1/sqrt(2)))./ corn_to_center_pix_dist; %this is mm/pixel



%% 2. Find velocities of every tag

Vel_cent = (Def_image_final_arranged(2:3,:) - Def_image_init_arranged(2:3,:))*(FPS/(F_final-F_init)); %Find velocities of the center of tags

%Now convert this from pixels/second to mm/second using tag_edge_length,
%Make this local to every tag to account for any small distortions

Vel_cent = Vel_cent .* repmat(pixel_mm_conv,2,1); %this repmat is needed as vel_cent has 2 comp that locally have same conversion as ith element of pixel...

%% 3. Find angular velocities of every tag

r_n_i = Def_image_init_arranged(4:end,:) - repmat(Def_image_init_arranged(2:3,:),4,1); %initial radius to corners for every tag

%convert to mm
r_n_i = r_n_i .* repmat(pixel_mm_conv,size(r_n_i,1),1);



r_n_f = Def_image_final_arranged(4:end,:) - repmat(Def_image_final_arranged(2:3,:),4,1); %final radius to corners for every tag

%convert to mm
r_n_f = r_n_f .* repmat(pixel_mm_conv,size(r_n_i,1),1);


s = (r_n_f - r_n_i)*(FPS/(F_final-F_init)); %velocity





omega_each = r_n_i([1 3 5 7],:).*s([2 4 6 8],:) - r_n_i([2 4 6 8],:).*s([1 3 5 7],:); %This is omega for each corner

%Now average omega:

omega_avg = mean(omega_each ,1);


%make frame vector
frame_vec = repmat(F_init,1,size(Def_image_final_arranged,2));

vel_prof = [Def_image_final_arranged(1,:);...
            Vel_cent;...
            omega_avg;...
            frame_vec];




end