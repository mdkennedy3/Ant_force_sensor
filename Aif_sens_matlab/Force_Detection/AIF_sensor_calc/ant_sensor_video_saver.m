%Function for recording video of deflection with data on it

function [str] = ant_sensor_video_saver(video, def_centers,Cal_img_info,Cal_img,Name_of_movie)



writerObj = VideoWriter([Name_of_movie,'.avi']);
open(writerObj);
%%%%%%%%%%%%%MOVIE%%%%%%%%%%%%%%%%%%%%%%%


h = figure('position',[120 120 1220 720]);


% for i = 2:nFrames+1
for i = 2:size(video,1)
    
    
    figure(h)
        
        subplot('position',[0,0.125,.5,.75])
        imshow(Cal_img)
        hold on
        
        
        
        
        tags = Cal_img_info{i}(1,:);
        p1 = Cal_img_info{i}(2:3,:);
        p2 = Cal_img_info{i}(4:5,:);
        p3 = Cal_img_info{i}(6:7,:);
        p4 = Cal_img_info{i}(8:9,:);
        p5 = Cal_img_info{i}(10:11,:);
        
        desire = size(Cal_img_info{i},2);
        for j = 1:desire
            
            plot(p1(1,j),p1(2,j),'p')
       
        str1 = {['T',num2str(Cal_img_info{i}(1,j))]};
        text(p1(1,j)+10,p1(2,j)+10,str1)
       
            plot(p2(1,j),p2(2,j),'o')
            plot(p3(1,j),p3(2,j),'o')
            plot(p4(1,j),p4(2,j),'o')
            plot(p5(1,j),p5(2,j),'o')
            
            
%             line([p1(1,j),p1(1,j)+def_centers{i}(2*(i-1)-1,j)],[p1(2,j),p1(2,j)+def_centers{i}(2*(i-1),j)],'Color',[1 0 0])
            line([p1(1,j),p1(1,j)+def_centers{i}(1,j)],[p1(2,j),p1(2,j)+def_centers{i}(2,j)],'Color',[1 0 0])
%             line([p1(1,j),p1(1,j)+sensor_def_centers_data{j+1,i}(1,1)],[p1(2,j),p1(2,j)+sensor_def_centers_data{j+1,i}(2,1)],'Color',[1 0 0])
        end
        title('Calibration')
        hold off
        
        
        %% Deformed image
        %     subplot(1,2,2)
        subplot('position',[.5,.125,.5,.75])
        imshow(video{i,1})
        hold on
        plot(video{i,2}(2,:), video{i,2}(3,:),  'p')
        
        %fix:
        for j=1:length(video{i,2}(2,:))
        str1 = {['T',num2str(video{i,2}(1,j))]};
        text(video{i,2}(2,j)+10,video{i,2}(3,j)+10,str1)
        end
        
        plot(video{i,2}(4,:), video{i,2}(5,:),  'o')
        plot(video{i,2}(6,:), video{i,2}(7,:),  'o')
        plot(video{i,2}(8,:), video{i,2}(9,:),  'o')
        plot(video{i,2}(10,:),video{i,2}(11,:), 'o')
        
        title('Deflection Video')
        hold off
    pause(.5)
%%%%%%%%%%%%%MOVIE%%%%%%%%%%%%%%%%%%%%%%%        
        frame = getframe(h);
        writeVideo(writerObj,frame);
%%%%%%%%%%%%%MOVIE%%%%%%%%%%%%%%%%%%%%%%%  
end

close(writerObj)

str = 'video recorded has been saved'

end