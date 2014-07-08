%% Sensor Deflection from Video
%Description:
%This code is for tracking deflection of sensors from video

%Developed by: Monroe Kennedy
%Date: 7/29/13
%Reference: April Tags code: http://people.csail.mit.edu/kaess/apriltags/


clear all
clc
close all

%This is for the undistortion code
% addpath(genpath('mex/'))


%% Place movie file, and calibration image here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath(genpath('~/Dropbox/Penn-ASU_share/'))

        Name_of_movie = '131009closeup';
        movie_extension = '.mp4';
        calibrated_image = '131009closeup_05.jpg';
        fc = [ 2668.28455   2691.37452 ];
        cc = [ 1442.37175   742.61228 ];
        
        Calib_matrix = eye(3,3).*repmat([fc,0]',1,3) + [zeros(3,1),zeros(3,1),[cc,1]']; %builds matrix with
        
        %Distortion:
        kc = [ -0.24121   0.08278   -0.00533   0.00062  0.00000 ];
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%12/30/13%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/Dropbox/Penn-ASU_share/Newformatvideo'))
Name_of_movie = 'video1080_60p';
movie_extension = '.mp4';
calibrated_image = 'M1010001.JPG';
fc = [ 1   1 ];
cc = [ 0   0 ];

Calib_matrix = eye(3,3).*repmat([fc,0]',1,3) + [zeros(3,1),zeros(3,1),[cc,1]']; %builds matrix with

%Distortion:
kc = [ 0  0   0   0  0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


obj = VideoReader([Name_of_movie,movie_extension])
Cal_img = imread(calibrated_image);



%% Setup:


tag_edge_length = 5; %Length of tags in mm



Cal_img = rgb2gray(Cal_img);

%Undistort Cal_img
Cal_img = undistort(Cal_img, Calib_matrix,  kc); %

[tags p1 p2 p3 p4 p5] = dapriltag(Cal_img);
Cal_info = [tags;...
    p1;...
    p2;...
    p3;...
    p4;...
    p5];




%% This is for extracting video frames

%Place the tags of interest here
%it is assumed there will be between 2 - 7 tags, with the last tag in the
%list being the center tag
%{
Tags_of_interest = [59366215950;...   %tag center
                    13868794921;...   %tag1
                    38943287133;...   %tag2
                    32013986588;...   %tag3
                    37856557532];%... %tag4
%}
Tags_of_interest = [2;...   %tag center
                    98;...  %tag1
                    62;...  %tag2
                    50;...  %tag3
                    38];    %tag4

%% Find the number of frames, and find acceptable images for comparison

dx_ends = [];
def_centers = cellmat;
Cal_img_info = cellmat;

%more accurate
%frame_count = read(obj); %This is apparently more accurate than .NumberofFrames, however takes longer     nFrames = obj.NumberOfFrames;
% nFrames = size(frame_count,4);

nFrames = obj.NumberOfFrames;
divise = 9;%1;%was 9 %35,7,5,51 => 1785
nFrames = floor(nFrames/divise);


video = cellmat;

sensor_data = cellmat;

%for velocity
Tags_of_interest_init = [];
flag = 0; %flag to indicate first velocity has been found

vel_prof = cellmat;
acc_prof = cellmat;

for i = 1:nFrames
    
    I = read(obj,floor(i*divise));
    I = rgb2gray(I);
    
    
    %% Undistort Image
    %undistort I and calimage
    %sensor.img = undistort(Image_grey, Cal_mat, distortion_coeff);
    I = undistort(I, Calib_matrix,  kc); %
    
    
    
    
    %%
    [tags p1 p2 p3 p4 p5] = dapriltag(I);
    
    %%%%%%%%%%%%%check to see if tags of interest are present%%%%%%%%%%%
    sensor_present = [];
    sensor = [];
    
    for k = 1:length(Tags_of_interest)
        if ~isempty(find(tags == Tags_of_interest(k)))
            sensor_present(end+1) = find(tags == Tags_of_interest(k));
            sensor(end+1) = tags(sensor_present(end));   %note that sensor is aligned as Tags_of_interest, but may not have all in set, but is ordered
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (~isempty(tags) && length(sensor_present) >= 2 && sensor(1) == Tags_of_interest(1))
        
        
        Def_info = [tags;...
            p1;...
            p2;...
            p3;...
            p4;...
            p5];
        
        %% Place velocity code here
        
        %Send in FPS in function
        FPS = 29.97;
        F_init = i-1;
        F_final = i;
        
        if flag == 0
            Tags_of_interest_init = sensor;
            def_img_info_init_frame = Def_info;
            flag = 1;
        else
            
            %update tags of interest
            Tags_of_interest_final = sensor;
            def_img_info_final_frame = Def_info;
            
            %Find velocity profile which include vel of cent tag, ang vel,
            %and corresponding frame
            [velocity_prof] = velocity_profile(def_img_info_init_frame,F_init, def_img_info_final_frame,F_final, tag_edge_length, Tags_of_interest_init,Tags_of_interest_final, FPS);
            
            if size(vel_prof,2) == 1 && size(vel_prof{1},1) == 1
                vel_prof{end} = velocity_prof;
            else
                vel_prof{end+1} = velocity_prof;
            end
            
            %now set init tags to latest
            Tags_of_interest_init = sensor;
            def_img_info_init_frame = Def_info;
            
        end
        
     
        %% Internal deflection
     
        video{end+1,1} = I;
        
        dx_ends = [];
        
        [dx_ends(end+1:end+2,:),def_centers{end+1}, Cal_img_info{end+1},Def_im_back] = deform_sens(Cal_img, Cal_info, Def_info, tag_edge_length, sensor);
        
        video{end,2} = Def_im_back;
        
        sensor_data{1,end+1} = i*divise;%This is the frame being looked at in the movie
        
        filler = 0; %this variable keeps track of how many times you filled a row
        for j = 2:size(sensor,2)
            
            index = find(Tags_of_interest == sensor(j)); %found where the first tag is (other than center in tags of interest)
            
            if sensor(j) ~= Tags_of_interest(j+filler)
                
                for k = j:(index-1)
                    sensor_data{j,end} = [];
                    filler = filler + 1;
                end
                
                %then add that value
                sensor_data{j,end} = dx_ends(:,j);
                
            else
                sensor_data{j,end} = dx_ends(:,j);
            end
        end
        
    end
    percent_complete = i/(nFrames)*100
end



%% Now find the Net Velocity profile

%Make matrix of velocity
velocity_prof_temp_frames = cell2mat(vel_prof);

tag_velocities_separate = cellmat;

for i = 1:length(Tags_of_interest)
    tag_velocities_separate{1,i} = Tags_of_interest(i);
end


for i = 1:size(velocity_prof_temp_frames,2)
    
    ind = find(Tags_of_interest' == velocity_prof_temp_frames(1,i)); %this determines which tag column this should belong to
    
    if size(tag_velocities_separate,1) == 1
        tag_velocities_separate{2,ind} = [velocity_prof_temp_frames(:,i)];
    else
        tag_velocities_separate{2,ind}(:,end+1) = [velocity_prof_temp_frames(:,i)];
    end
end

%% Now give graph for deflection over frames
[str] = plot_sensor_data(sensor_data,tag_velocities_separate,divise,Name_of_movie,'median_filt') %Type is type of filtering: 'impulse' or 'median_filt'

%% Save video:
[str] = ant_sensor_video_saver(video, def_centers,Cal_img_info,Cal_img,Name_of_movie)

%% Save data if desired:
save([Name_of_movie, '_sensor_data'],'sensor_data','-v7.3')
save([Name_of_movie, '_sensor_vel_data'],'tag_velocities_separate','-v7.3')
