function [str] = plot_sensor_data(sensor_data,tag_velocities_separate,divise, Name_of_movie,type)

if strcmpi(type,'median_filt')
   
    close all
    for k = 2:size(sensor_data,1)
        figure(k)
        
        hold on
        sens_data_temp = [sensor_data{k,:}];
        
        sens_data_frame = [];
        
        %Find essential frames for comparison
        for j = 1:size(sensor_data,2)
            if ~isempty(sensor_data{k,j})
                sens_data_frame(end+1) = sensor_data{1,j};
            end
        end

        %Set window size: 
        %250
        win_size = 10; %this means use median of 10 data points
        
        xwindow = zeros(1,length(sens_data_temp(1,:)));
        ywindow = zeros(1,length(sens_data_temp(2,:)));
        
        for i = 1:length(sens_data_temp(1,:))
        
            if i <= win_size
        xwindow(i) = sens_data_temp(1,i);
        ywindow(i) = sens_data_temp(2,i);
            else
        xwindow(i) = median(sens_data_temp(1,(i-win_size):i));
        ywindow(i) = median(sens_data_temp(2,(i-win_size):i));
                
            end
        end

        
%KEEP BELOW:
        plot(sens_data_frame(:),sens_data_temp(1,:),'k',sens_data_frame(:),xwindow,'-.',...
            sens_data_frame(:),sens_data_temp(2,:),'r',sens_data_frame(:),ywindow,'--')
         legend('x position','x position median filtered','y position','y position median filtered')        
       
        title(['Position of Leg ',num2str(k-1)])
        xlabel('frames')
        ylabel('mm')
        
        hold off
        
    h = gcf    
    filename = [Name_of_movie,'_Sensor_leg',num2str(k-1)]
     print(h,'-dpdf', '-r200',filename) 
        
    end
    
     save([Name_of_movie, '_sensor_data'],'sensor_data','-v7.3')
    
    str = 'Data saved'
   
    
    
    
    
elseif strcmpi(type, 'impulse')
    
    
    close all
    window_size = 10;
    fir = ones(1,window_size);
    
    for k = 2:size(sensor_data,1)
        figure(k)
        
        hold on
        sens_data_temp = [sensor_data{k,:}];
        
        
        sens_data_frame = [];
 
        %Find essential frames for comparison
        for j = 1:size(sensor_data,2)
            if ~isempty(sensor_data{k,j})
                sens_data_frame(end+1) = sensor_data{1,j};
            end
        end
        
        
        %Filter        
        xwindow = conv(sens_data_temp(1,:), fir/window_size); xwindow = xwindow(1:length(sens_data_frame(:)));
        ywindow = conv(sens_data_temp(2,:), fir/window_size); ywindow = ywindow(1:length(sens_data_frame(:)));
        
        
        %Plot just filter
        plot(sens_data_frame(:),xwindow,'-.',...
             sens_data_frame(:),ywindow,'--')
        legend('x position impulse filtered','y position impulse filtered')
        
        %Plot both filter and data
%         plot(sens_data_frame(:),sens_data_temp(1,:),'k',sens_data_frame(:),xwindow,'-.',...
%             sens_data_frame(:),sens_data_temp(2,:),'r',sens_data_frame(:),ywindow,'--')
%          legend('x position','x position impulse filtered','y position','y position impulse filtered')
        
        
        xlabel('frames')
        ylabel('mm')
        title(['Position of Leg ',num2str(k-1)])
        hold off
    end
    
    
    
    save([Name_of_movie, '_sensor_data'],'sensor_data','-v7.3')
    
    str = 'Data saved'
    
else    %Right now default to impulse filter
    
    close all
    window_size = 10;
    fir = ones(1,window_size);
    
    for k = 2:size(sensor_data,1)
        figure(k)
        
        hold on
        sens_data_temp = [sensor_data{k,:}];
        
        
        sens_data_frame = [];
        
        for j = 1:size(sensor_data,2)
            if ~isempty(sensor_data{k,j})
                sens_data_frame(end+1) = sensor_data{1,j};
            end
        end
        
        
        
        xwindow = conv(sens_data_temp(1,:), fir/window_size); xwindow = xwindow(1:length(sens_data_frame(:)));
        ywindow = conv(sens_data_temp(2,:), fir/window_size); ywindow = ywindow(1:length(sens_data_frame(:)));
        
        
        plot(sens_data_frame(:),sens_data_temp(1,:),'k',sens_data_frame(:),xwindow,'-.',...
            sens_data_frame(:),sens_data_temp(2,:),'r',sens_data_frame(:),ywindow,'--')
        legend('x position','x position filtered','y position','y position filtered')
        xlabel('frames')
        ylabel('mm')

        title(['Position of Leg ',num2str(k-1)])
        hold off
    end
    
    
    
    save([Name_of_movie, '_sensor_data'],'sensor_data','-v7.3')
    
    str = 'Data saved'
    
end


%% Plot tag velocities vs. Frame

%make plots for each tag
for k = 1:size(tag_velocities_separate,2)
    
    frame_ind = tag_velocities_separate{2,k}(end,:).*divise; %frame indicies
    figure
    x = plot(frame_ind, tag_velocities_separate{2,k}(2,:),'r',frame_ind, tag_velocities_separate{2,k}(3,:),'b')
    legend([x],'X velocity','Y velocity')
    title(['Velocity of tag ',num2str(tag_velocities_separate{2,k}(1,1))])
    xlabel('Frame')
    ylabel('velocity: mm/s')
    
    h = gcf    
    
    if k == 1
    filename = [Name_of_movie,'_Vel_Sensor_center']
    else
    filename = [Name_of_movie,'_Vel_Sensor_leg',num2str(k-1)]
    end
    print(h,'-dpdf', '-r200',filename) 
    
end




end