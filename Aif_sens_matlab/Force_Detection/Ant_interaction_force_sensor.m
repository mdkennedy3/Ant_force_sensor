%Ant Interaction Force Sensor (AIFS)
%UPENN & ASU

clear all
clc
close all


addpath(genpath('AIF_sensor_calc/'))


%% Camera Calibration

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


%% Get sensor data for experiment

%%%%%%%Setup folders with data of interest%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%%%%%%%%%%%%%%%%12/30/13%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('~/Dropbox/Penn-ASU_share/Newformatvideo'))
Name_of_movie = 'video1080_60p';
movie_extension = '.mp4';
calibrated_image = 'M1010001.JPG';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%{
%%%%%%%%%%%%%%%%%test_for_center_trajectory%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('~/Dropbox/Penn-ASU_share/131009_closeup'))
Name_of_movie = '131009closeup';
movie_extension = '.mp4';
calibrated_image = '131009closeup_03.jpg';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%Testing force calibration routines
%{
%%%%%%%%%%%%%%%%5/30/14%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('~/Dropbox/Penn-ASU_share/140506/')
Name_of_movie = 'verifT1';
movie_extension = '.mp4';
calibrated_image = 'Calibration_image.JPG';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%%%%%%%%%%%%%%%%%7/8/14%%%%%%%%%%%%%%%%%%%%%%%%%%%
 addpath('~/Dropbox/Penn-ASU_share/140612/') %
 Name_of_movie = '00007';
 movie_extension = '.mp4';
 calibrated_image = 'S2400053.JPG';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tag_edge_length = 5; %Length of tags in mm
Frame_rate_collect = 9; %This basically determines how many frames to skip between collecting data
FPS = 25; %Frame rate of camera

plot_and_save_tag_data = 0; %plot and save position and velocity data
plot_cent_trajectory = 0; %plot the trajectory of the center tag
save_tag_pos_video = 0; %save video showing deflection on calibrated image next to deflection video
save_tag_pos_data = 0; %save tag position data
save_tag_vel_data = 0; %save tag velocity data
save_cent_trajectory = 0; %save the trajectory of the center
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pos_prof, vel_prof, cent_traj] = aif_sensor(Name_of_movie, movie_extension, calibrated_image, Calib_matrix, kc,...
    FPS, tag_edge_length, Frame_rate_collect, plot_and_save_tag_data,...
    save_tag_pos_video, save_tag_pos_data, save_tag_vel_data,save_cent_trajectory,plot_cent_trajectory);

%note for cent_traj that pixels in matlab have (0,0) in top left corner of frame, hence y will be 'flipped' when plotted
%% Get forces from deflections


%%%%%%%%%% First Load constant matricies for all tags (saved from
%%%%%%%%%% calibration routine)
disp('Make Sure you put conts matricies in this folder or addpath to access them')


Fmag_tags = cellmat
Fx_tags = cellmat
Fy_tags = cellmat

%% Find Forces for Tag1  (Comment out if you don't have)
 Tag_of_interest = 1;
 %Load what exists
load const_mat_T1.mat
[Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof, Cont_mat, Tag_of_interest);
Fmag_tags{1} = Fmag;
Fx_tags{1} = Fx;
Fy_tags{1} = Fy;

%%%%%%%%%%%%%%%%%Plot and save if desired%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%% I've commented out the below as for this test script only Tag 1 is present
%{
%% Find Forces for Tag2 (Comment out if you don't have)
Tag_of_interest = 2;
Load what exists
load const_mat_T2.mat
[Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof, const_mat_T2, Tag_of_interest);
Fmag_tags{2} = Fmag;
Fx_tags{2} = Fx;
Fy_tags{2} = Fy;
%%%%%%%%%%%%%%%%%Plot and save if desired%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

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

%% Find Forces for Tag3 (Comment out if you don't have)
Tag_of_interest = 3
%Load what exists
load const_mat_T3.mat
[Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof, const_mat_T3, Tag_of_interest);
Fmag_tags{3} = Fmag;
Fx_tags{3} = Fx;
Fy_tags{3} = Fy;
%%%%%%%%%%%%%%%%%Plot and save if desired%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

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
%% Find Forces for Tag4 (Comment out if you don't have)
Tag_of_interest = 4;
%Load what exists
 load const_mat_T4.mat
[Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof, const_mat_T4, Tag_of_interest);
Fmag_tags{4} = Fmag;
Fx_tags{4} = Fx;
Fy_tags{4} = Fy;
%%%%%%%%%%%%%%%%%Plot and save if desired%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

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


%}

%% Now Save force Data
save(['Fmag_tags'],'Fmag_tags','-v7.3');
save(['Fx_tags'],'Fx_tags','-v7.3');
save(['Fy_tags'],'Fy_tags','-v7.3');



