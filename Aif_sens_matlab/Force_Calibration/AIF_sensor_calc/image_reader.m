% This is testing the april code file
close all
clc
clear

%addpath('deflectionpics/')
addpath('test_images/')
%I = imread('photo.JPG');
%I = imread('find1.jpeg');
%I = imread('find2.JPG');
addpath('~/Dropbox/Penn-ASU_share/140506/')
%I = imread('tags.JPG');

%% Import Image
M = cellmat;

M{1} = imread('Calibration_image.JPG')


%{
M{1} = imread('131009closeup_01.jpg')
M{2} = imread('131009closeup_02.jpg')
M{3} = imread('131009closeup_03.jpg')
M{4} = imread('131009closeup_04.jpg')
M{5} = imread('131009closeup_05.jpg')
M{6} = imread('131009closeup_06.jpg')
M{7} = imread('131009closeup_07.jpg')
M{8} = imread('131009closeup_08.jpg')
M{9} = imread('131009closeup_09.jpg')
M{10} = imread('131009closeup_10.jpg')
M{11} = imread('131009closeup_11.jpg')
%}
%%

for i = 1:size(M,2)

I = M{i}; 

I = rgb2gray(I);

[tags p1 p2 p3 p4 p5] = dapriltag(I);


information = [tags;...
                p1;...
                p2;...
                p3;...
                p4;...
                p5];
      
%Now draw the corners on the screen



figure(i)
imshow(I)
hold on


desire = size(information,2);

for i = 1:desire

 plot(p1(1,i),p1(2,i),'p')
 plot(p2(1,i),p2(2,i),'o')
 plot(p3(1,i),p3(2,i),'o')
 plot(p4(1,i),p4(2,i),'o')
 plot(p5(1,i),p5(2,i),'o')

 %average pt
cent = mean([p2(1:2,i),p3(1:2,i),p4(1:2,i),p5(1:2,i)],2);

plot(cent(1),cent(2),'r')
end


%annotation('arrow',[.39 .39],[.11 .39])

hold off

information(:,1:desire);



end
%% Tag information

%tag 1: 13868794921
%tag 2: 38943287133
%tag 3: 32013986588
%tag 4: 37856557532
%tag 5: 47628935088
%tag 6: 7661251159
%tag center: 59366215950

%% This is for extracting video frames
% 
% %obj = mmreader('tag_movie.MOV');
% %video = obj.read();
% 
% obj = VideoReader('tag_movie.mp4');
% 
% nFrames = obj.NumberOfFrames;
% 
% %video = obj.read();
% 
% video = cellmat;
% 
% for i = 1:nFrames
% 
%     I = read(obj,i);
%     I = rgb2gray(I);
%  
%     [tags p1 p2 p3 p4] = dapriltags(I);
% 
% information = [tags;...
%                 p1;...
%                 p2;...
%                 p3;...
%                 p4];
%    
%             video{end+1,1} = I;
%             video{end,2} = information;
%             
%             
% % figure(1)            
% % imshow(I)
% % hold on
% % 
% % plot(p1(1,:),p1(2,:),'p')
% % plot(p2(1,:),p2(2,:),'o')
% % plot(p3(1,:),p3(2,:),'o')
% % plot(p4(1,:),p4(2,:),'o')
% % 
% % hold off
%             
%             
%             
% end
%    
% 
% 
% 
% for i = 2:nFrames+1
%    
%     figure% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 



%     imshow(video{i,1})
%     hold on
%      
%     plot(video{i,2}(2,:),video{i,2}(3,:),'p')
%     plot(video{i,2}(4,:),video{i,2}(5,:),'p')
%     plot(video{i,2}(6,:),video{i,2}(7,:),'p')
%     plot(video{i,2}(8,:),video{i,2}(9,:),'p')
% end
% 
% 
% 
% 
% 
















