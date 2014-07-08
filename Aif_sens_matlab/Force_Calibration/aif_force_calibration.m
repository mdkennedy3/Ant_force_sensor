%This Script is for force calibration, this script saves vector matricies
%that capture Nonlinear relationship between force and deflection

%This script must be run with new calibration videos for each tag as the
%focus of the video, and indicated in the second portion when prompted

clear all
clc
close all

%Needed to access control scripts
addpath(genpath('AIF_sensor_calc/'))



%% Camera Calibration

%Add appropriate path for calibration matrix, and/or load it below

%Note: with Matlab 2014 higher, there is a built in >> cameraCalibrator
% using this given images of checkerboard as in http://www.vision.caltech.edu/bouguetj/calib_doc/
% For the Camera Calibrator software you must specify the edge length
%of checkerboard square (e.g. 12mm) then let it calculate.
% the camera calibrator software will produce calibration matrix (3x3) and
% distortion matrix which is composed of two: radial (3 components) and
% tangential (2 components)
% The distortion coefficient must be kc = [rad(2x), tang(2x), 0];
% The calibration matrix will be called 'intrinsic_matrix'

%Step 1: place calibration images in one folder: e.g. Cam_Calibration
%Step 2: run >> cameraCalibrator
%Step 3: click on 'Add Images', navigate to folder (e.g. Cam_Calibration)
%Step 4: make sure radial distortion coefficients is on 2 coefficients
%Step 5: click "Calibrate"
%Step 6: click "Export Camera Parameters"
%Step 7:

%adjust accordingly:
%load cam_param_131009_c.mat

camera_calib_2014 = exist('cameraParameters','var');

if camera_calib_2014 == 1
    Calib_matrix = cameraParameters.IntrinsicMatrix'; %transpose as last column should not contain 0's
    kc = [cameraParameters.RadialDistortion, cameraParameters.TangentialDistortion, 0];
else
    %default
    fc = [ 1   1 ];
    cc = [ 0   0 ];
    Calib_matrix = eye(3,3).*repmat([fc,0]',1,3) + [zeros(3,1),zeros(3,1),[cc,1]']; %builds matrix with
    %Distortion:
    kc = [ 0  0   0   0  0 ];
end



%% Sensor Calibration

%%%%%%%Setup folders with data of interest%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%5/30/14%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('~/Dropbox/Penn-ASU_share/140506/')
% Name_of_movie = 'verifT1';
% movie_extension = '.mp4';
% calibrated_image = 'Calibration_image.JPG';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%5/30/14%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('~/Dropbox/Penn-ASU_share/140612/')
Name_of_movie = '00009';
movie_extension = '.mp4';
calibrated_image = 'S1.JPG';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tag_edge_length = 5; %Length of tags in mm
Frame_rate_collect = 5; %This basically determines how many frames to skip between collecting data
FPS = 25; %Frame rate of camera


%Data Saving Options
plot_and_save_tag_data = 0; %plot and save position and velocity data
save_tag_pos_video = 1; %save video showing deflection on calibrated image next to deflection video
save_tag_pos_data = 0; %save tag position data
save_tag_vel_data = 0; %save tag velocity data
save_cent_trajectory = 0;
plot_cent_trajectory =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[pos_prof_sens_calibrate, vel_prof_sens_calibrate, cent_traj] = aif_sensor(Name_of_movie, movie_extension, calibrated_image, Calib_matrix, kc,...
    FPS, tag_edge_length, Frame_rate_collect, plot_and_save_tag_data,...
    save_tag_pos_video, save_tag_pos_data, save_tag_vel_data,save_cent_trajectory,plot_cent_trajectory);


%% Force Constant Matrix

%If desired comment out below prompt section and enter data manually into
%masses and index_of_masses variables

%% Automation Section
disp('At this point look at video and determine times that weights were changed, and with margin, ')
disp(' figure out the time sections that coorespond to the weights.')
disp('    ')
disp(['The data has approximately ', num2str(size(pos_prof_sens_calibrate,2)) ,' points, hence given the times and FPS: ', num2str(FPS), ' this script will determine the appropriate indicies for each mass'])
%This is approximate, because tag may not be present in all
disp('    ')
disp('It is wise to find the times with a little buffer to avoid overlap in calibration, for instance if video has 1-5sec mass 1, 5-10sec mass2, then say 1 -4.8sec is mass 1, and 5.2-10 sec is mass 2.')
disp('    ')
disp('Either place the time in prompts below, or manually change code (do calculation to find index) and skip prompts. ')
disp('    ')
prompt = 'How many weights sections are in the video? (enter scalar e.g. 4): ';
str_num_of_weights = input(prompt,'s')
disp('    ')
num_of_weights = str2num(str_num_of_weights);



disp('   ')

prompt ='Please indicate the tag of interest: 1 for tag 1, 2 for tag 2, etc: ';
disp('   ')

tag_of_interest_str = input(prompt,'s');

Tag_of_interest = str2num(tag_of_interest_str); %Save the indicated tag of interest to str



masses = []; %mass in grams
%e.g. masses = [ 1.57, 1.16, 0.97,  0.75 ]; %grams
for i = 1:num_of_weights
    prompt= [' Enter the weight in grams of the ', num2str(i) ,' Mass: '  ];
    str_latest_mass = input(prompt,'s')
    latest_mass = str2num(str_latest_mass);
    masses(1,end+1) = latest_mass;
end
disp(' Thank you, ')
disp('    ')
disp(' Now we will enter the indexes for each mass: ')
disp('    ')

index_of_masses = []; %indexes of masses e.g.:
% index_of_new_mass = [1, 122;...  %indexes for first section
%                      153, 274;...%indexes for second section
%                      305, 457;...%indexes for third section
%                      488, size(Tags_bound_cond_sens_calibrate,1)]; %indexes for fourth section

for i = 1:num_of_weights
    if i == 1
       prompt = 'For first weight please indicate time (secs) when it stops in calibration video created (give small buffer: 5.25sec, enter 5sec), for instance mass_1 time 0:5 sec,  so enter 5: ';
        time_index_first_str =  input(prompt,'s')
        time_index_first = str2num(time_index_first_str);
        index_of_masses(end+1,:) = [1, floor(FPS*time_index_first)];
        
    elseif i == num_of_weights
        prompt = 'For last weight please indicate time(sec) when it begins, for instance m1: time 16.5 : 20 sec), so enter 16.5: ';
        time_index_last_str =  input(prompt,'s');
        time_index_last = str2num(time_index_last_str);
        
        %ensure that the last index you input for the indexes actually
        %exist: for instance the end of pos_prof may only contain info on 1
        %tag, and so need to find last index that contains tag info
        find_last = [];
        count = size(pos_prof_sens_calibrate,2);
        while isempty(find_last) && (count ~= 1)
            if ~isempty(pos_prof_sens_calibrate{Tag_of_interest+1, count})
                find_last = count;
            end
            count = count - 1;
        end
        
        
        index_of_masses(end+1,:) = [floor(FPS*time_index_last), find_last];
    else
        prompt = ['Please indicate the 1st time (sec) of the ', num2str(i), ' mass. E.g. time 5.25 : 8.5 sec, enter now the scalar 5.25: '];
        disp('    ')
        time_first_index_ith_str  = input(prompt,'s')
        time_first_index_ith = str2num(time_first_index_ith_str);
        prompt = ['Please indicate the 2nd time (sec) of the ', num2str(i), ' mass. E.g. time 5.25 : 8.5 sec, enter now the scalar 8.5: '];
        disp('    ')
        time_second_index_ith_str  = input(prompt,'s')
        time_second_index_ith = str2num(time_second_index_ith_str);
        
        index_of_masses(end+1,:) = [floor(FPS*time_first_index_ith), floor(FPS*time_second_index_ith)];
    end
end

disp('   ')

disp('Here are masses and indexes stored:')
masses
index_of_masses



%% %%%%%%% Manual Input for masses/index_of_masses %%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Set Masses hung in video
%masses = [ 1.57, 1.16, 0.97,  0.75 ]; %grams

%%%%%% Save Video above and
%index_of_masses = [1, 122;...  %indexes for first section
%                    153, 274;...%indexes for second section
%                   305, 457;...%indexes for third section
%                  488, size(Tags_bound_cond_sens_calibrate,1)]; %indexes for fourth section

%Manually Set Tag of Interest
% Tag_of_interest = 1; %specify tag of interest

%% Now Caculate the Constant Matrix vector and save

[Cont_mat] = aif_elliptic_force_surf(masses, pos_prof_sens_calibrate,index_of_masses , Tag_of_interest);


%now save const_matrix
if Tag_of_interest == 1
save(['const_mat_T1'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

elseif Tag_of_interst == 2

save(['const_mat_T2'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

elseif Tag_of_interest == 3

save(['const_mat_T3'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

elseif Tag_of_interest == 4

save(['const_mat_T4'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

elseif Tag_of_interest == 5
save(['const_mat_T5'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

elseif Tag_of_interest == 6
save(['const_mat_T6'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)

else
save(['const_mat'],'Cont_mat','-v7.3'); %This is save('name_to_be_saved_as','current_variable','-v7.3)
    
end



%% Now immediately test and see:
%{
[Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof_sens_calibrate, Cont_mat, Tag_of_interest);

index_count = 1:length(Fmag);
figure
plot(index_count,Fmag)
title('Force Magnitude')
xlabel('Data Index')
ylabel('Force mN')
h = gcf
print(h,'-dpdf', '-r200','Fmag_6_7')

figure
plot(index_count,Fx)
h = gcf
title('Fx Magnitude')
xlabel('Data Index')
ylabel('Force mN')
print(h,'-dpdf', '-r200','Fx_6_7')

figure
plot(index_count,Fy)
title('Fy Magnitude')
xlabel('Data Index')
ylabel('Force mN')
h = gcf
print(h,'-dpdf', '-r200','Fy_6_7')

figure
plot(index_count,Fy,'r',index_count,Fx,'b')
title('Fx and Fy Magnitudes')
xlabel('Data Index')
ylabel('Force mN')
legend('Fy','Fx')
h = gcf
print(h,'-dpdf', '-r200','Fy_and_Fx_6_7')
%}
