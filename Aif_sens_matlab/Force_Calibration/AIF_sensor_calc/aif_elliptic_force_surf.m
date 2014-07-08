function [Const] = aif_elliptic_force_surf(masses, pos_prof, index_data, Tag_of_interest)
%This function used for force calibration takes in masses, and tag
%positions, and using the index trains a set of constants that fit a 3D
%ellipcital surface to the data, where the gradient is the direction of the
%force, and the 'height' is the magnitude of the force. We believe that
%this elliptical surface will appropriately capture the nonlinearity of the
%helical spring
%{
masses = [ 1.57, 1.16, 0.97,  0.75 ]; %grams

index_data  = [1, 122;...  %indexes for first section
               153, 274;...%indexes for second section
               305, 457;...%indexes for third section
               488, size(Tags_bound_cond_sens_calibrate,1)]; %indexes for fourth section
 
%}


%% Setup

%1. Find the force magnitudes in each section:
%F_mag = masses*1e-3*9.8; %convert mass from grams to kg, then to N
F_mag = masses*9.8; %convert mass g to mN



%2. make data matrix that repmats Forces and places them next to
%correponding deflections: D_i = [F_i, x_i, y_i]
Data = [];
for i = 1:size(index_data,1)
    for j = index_data(i,1):index_data(i,2)
        
        if ~isempty(pos_prof{1+Tag_of_interest,j})
           
        Data(end+1,:) = [F_mag(i), pos_prof{1+Tag_of_interest,j}(1),pos_prof{1+Tag_of_interest,j}(2) ]; % [mN, mm,mm]
        end
    end 
end


%3. Make calibration force vector matrix (using boundary conditions, to be
%used later to confirm gradient direction approach, when testing on
%training set)


%% Build Force Vector: [F_i]_{Nx1}
%Just take D_i(:,1)

F_vec = Data(:,1); %mN

%% Build Training Matrix: [x_i^2, x_i, 1, y_i^2, y_i, 1]_{Nx6}

ellip_mat = [];

for i = 1:size(Data,1)
    x_i = Data(i,2); %mm
    y_i = Data(i,3); %mm
ellip_mat(end+1,:) = [x_i^2, x_i, 1, y_i^2, y_i, 1];    
end


%% Invert to find Constants matrix

Const = pinv(ellip_mat)*F_vec;

end
