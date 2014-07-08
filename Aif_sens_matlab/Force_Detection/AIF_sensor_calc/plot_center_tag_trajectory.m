function [str] = plot_center_tag_trajectory(cent_pix_traj, movie_frame_height, movie_frame_width, tag_length,Name_of_movie)

%since center_pixel_traj contains [p1;p2;p3;p4;p5] for the center tag, then
% for every instant the correlation to mm can be found at every point based
% on p2,p3,p4,p5




    p2 = cent_pix_traj(3:4,1); %pos of corner 1
    p3 = cent_pix_traj(5:6,1); %pos of corner 2
    p4 = cent_pix_traj(7:8,1); %pos of corner 3
    p5 = cent_pix_traj(9:10,1); %pos of corner 4
     
    scale = tag_length/(mean([norm(p2-p3),norm(p3-p4), norm(p4-p5),norm(p5-p2)]));  %scale = mm/pixels
    

center_tag_position = [];
for i = 1:size(cent_pix_traj,2)
  
    %since in a movie in matlab 0,0 is in the top left corner, I must subract
    %the height from teh y position so that it is accurate when compared to the
    %movie
    
    p1 = [cent_pix_traj(1,i);movie_frame_height-cent_pix_traj(2,i)]; %position of center of tag
        
    center_tag_position = [center_tag_position,p1.*scale];
    
end


%not plot this and save plot

figure


axis([0,movie_frame_width*scale, 0, movie_frame_height*scale])
hold on


plot(center_tag_position(1,:),center_tag_position(2,:));
xlabel('center tag x position (mm)')
ylabel('center tag y position (mm)')
title('Trajectory of center tag')
hold off

h = gcf    
filename = [Name_of_movie,'_trajectory_of_center_tag']
print(h,'-dpdf', '-r200',filename) 


str = 'Data saved';


end