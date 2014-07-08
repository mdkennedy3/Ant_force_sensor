function [Consts] = aif_force_calibration(Tags_bound_cond, masses, index_of_new_mass, ws, Tag_of_interest)

%% input
%{
Tags_bound_cond:  Tags_bound_cond{:,1} = tag_ang_BC;
                  Tags_bound_cond{:,2} = Ls_mat; %This value will remain constant
                  Tags_bound_cond{:,3} = tag_init_position; %This value will remain constant
                  Tags_bound_cond{:,4} = tag_pos_bound_cond;

masses: [m1; m2; ... ; mN]; N masses in grams

index_of_new_mass: if there are N data points (N rows of Tags_bound_cond),
then the index where a new mass was hung (assuming video has been edited
such that for a specific tag

Note that for index_of_new_mass: for
index_of_new_mass = [1, 122;...  %indexes for first section  Hanging Mass 1
                     153, 274;...%indexes for second section Hanging Mass 2
                     305, 457;...%indexes for third section  Hanging Mass 3
                     488, 640]; %indexes for fourth section  Hanging Mass 4

ws: %weight/length of wire: mN/mm

Tag_of_interest: Specifiy the Tag that code should focus on for analysis.
Constants given will be specific for this tag.

%}
%% Setup
%First determine how many mass segments need to be analyzed and split the
%code accordingly


if (size(index_of_new_mass,1) == 1)  %hence only one mass was given
    
    F_mag = masses*1e-3*9.81; %convert grams mass to weight (Newtons)
    
    coil_weight_den = ws*1e3; %since ws is mN/mm and L is mm, then ws*L*1e3 is Newton weight of spring, we will not multiply by length here as its accounted for in Bern, Eul. eqns
    %Now it is assumed that the weight is pulling down in the frame.
    
    %% Make force vector
    Force_vec = []; %Force in Newtons, format: F(i,:) = [Fx, Fy]
    
    for i = 1:size(Tags_bound_cond,2)
        
        Fx_cal = -F_mag*sin(-Tags_bound_cond{i,1}(1));
        Fy_cal = -F_mag*cos(-Tags_bound_cond{i,1}(1));
        Force_vec(end+1,:) = [Fx_cal,Fy_cal]; %force from weight in calibrated frame
        
    end
    
    coil_weight = [];
    
    for i = 1:size(Tags_bound_cond,2)
        
        cw_x_cal = -coil_weight_den*sin(-Tags_bound_cond{i,1}(1));
        cw_y_cal = -coil_weight_den*cos(-Tags_bound_cond{i,1}(1));
        coil_weight(end+1,:) = [cw_x_cal,cw_y_cal]; %force from weight in calibrated frame
        
    end
    
    
    
    
    %% Build B.C. Matrix for N points
    bc_matrix = [];
    for i = 1:size(Tags_bound_cond,2)
        bc_matrix(end+1:end+4,1) = [0;... dy/dx(0), angle of center tag defined to be 0
            Tags_bound_cond{i,1}(1+Tag_of_interest);...  dy/dx(L), angle of Tag w.r.t. center tag (in center tag frame)
            0;... y(o), position of center tag, defined to be 0
            Tags_bound_cond{i,4}(Tag_of_interest)]; %y(L), position of Tag of interest
        
    end
    
    
    
    %% Build Var_matrix for N points
    Var_matrix = [];
    
    for i = 1:size(Force_vec,1)
        
        Fx = Force_vec(i,1);
        Fy = Force_vec(i,2);
        wx = coil_weight(i,1);
        wy = coil_weight(i,2);
        Ls = Tags_bound_cond{i,2}(Tag_of_interest);
        
%         
%         Var_matrix(end+1:end+4,1) = [ 0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  1, 0;...
%                                       0, 0, 0,           (wy*(Fx + Ls*wx)^3)/6, 0, 0,      ((Fx + Ls*wx)^2*(Fy + Ls*wy))/2, 0, 0,   (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx))/2,          0, 0,  (Ls^2*(3*Fy + Ls*wy))/6,  1, 0;...
%                                       0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  0, 1;...
%                  (wy*(Fx + Ls*wx)^4)/24, 0, 0, ((Fx + Ls*wx)^3*(Fy + Ls*wy))/6, 0, 0, (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx)^2)/4, 0, 0, (Ls^2*(3*Fy + Ls*wy)*(Fx + Ls*wx))/6, Fx + Ls*wx, 0, (Ls^3*(4*Fy + Ls*wy))/24, Ls, 1];
    
    
        Var_matrix(end+1:end+4,1) = [ 0,                                0,                                     0,                                     0,          0,                         0,  1, 0;...
                                      0,            (wy*(Fx + Ls*wx)^3)/6,       ((Fx + Ls*wx)^2*(Fy + Ls*wy))/2,    (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx))/2,          0,   (Ls^2*(3*Fy + Ls*wy))/6,  1, 0;...
                                      0,                                0,                                     0,                                     0,          0,                         0,  0, 1;...
                 (wy*(Fx + Ls*wx)^4)/24,  ((Fx + Ls*wx)^3*(Fy + Ls*wy))/6,  (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx)^2)/4,  (Ls^2*(3*Fy + Ls*wy)*(Fx + Ls*wx))/6, Fx + Ls*wx,  (Ls^3*(4*Fy + Ls*wy))/24, Ls, 1];
    
    
    end
    
    %% Solve for Consts Matrix
    
    Consts = pinv(Var_matrix)*bc_matrix;
    
    
    
else %multiple masses are given in this dataset, best case as the more wieghts are given the better trained the vector consts will be.
    
    F_mag = masses*1e-3*9.81; %convert grams mass to weight (Newtons)
    coil_weight_den = ws*1e-3; %since ws is mN/mm and L is mm, then ws*L*1e-3 is Newton weight of spring, we will not multiply by length here as its accounted for in Bern, Eul. eqns
    %Now it is assumed that the weight is pulling down in the frame.
    
    
    %% Make force vector, where Forces match segments where they hang, (using index_of_new_mass)
    
    Force_vec = []; %Force in Newtons, format: F(i,:) = [Fx, Fy]
%     count_position = 1;
    for j = 1:size(index_of_new_mass,1)
        for i = index_of_new_mass(j,1):index_of_new_mass(j,2)
            Fx_cal = -F_mag(j)*sin(-Tags_bound_cond{i,1}(1));
            Fy_cal = -F_mag(j)*cos(-Tags_bound_cond{i,1}(1));
            Force_vec(end+1,:) = [Fx_cal,Fy_cal]; %force from weight in calibrated frame
            
        end
%         count_position = index_of_new_mass(j);
    end
    
    coil_weight = [];
    count_position = 1;
    for j = 1:size(index_of_new_mass,1)
        for i = index_of_new_mass(j,1):index_of_new_mass(j,2)
            
            cw_x_cal = -coil_weight_den*sin(-Tags_bound_cond{i,1}(1));
            cw_y_cal = -coil_weight_den*cos(-Tags_bound_cond{i,1}(1));
            coil_weight(end+1,:) = [cw_x_cal,cw_y_cal]; %force from weight in calibrated frame
            
        end
        count_position = index_of_new_mass(j);
    end
    
    %% Build B.C. Matrix for N points
    bc_matrix = [];
    
    %HERE
%     for i = 1:size(Tags_bound_cond,1)
 for j = 1:size(index_of_new_mass,1)
        for i = index_of_new_mass(j,1):index_of_new_mass(j,2)
            
            bc_matrix(end+1:end+4,1) = [0;... dy/dx(0), angle of center tag defined to be 0
            Tags_bound_cond{i,1}(1+Tag_of_interest);...  dy/dx(L), angle of Tag w.r.t. center tag (in center tag frame)
            0;... y(o), position of center tag, defined to be 0
            Tags_bound_cond{i,4}(Tag_of_interest)]; %y(L), position of Tag of interest
        
        end
 end
    
    
    %% Build Var_matrix for N points
    Var_matrix = [];
    
    for i = 1:size(Force_vec,1)
        
        Fx = Force_vec(i,1);
        Fy = Force_vec(i,2);
        wx = coil_weight(i,1);
        wy = coil_weight(i,2);
        Ls = Tags_bound_cond{i,2}(Tag_of_interest);
        
        
        Var_matrix(end+1:end+4,:) = [ 0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  1, 0;...
            0, 0, 0,           (wy*(Fx + Ls*wx)^3)/6, 0, 0,      ((Fx + Ls*wx)^2*(Fy + Ls*wy))/2, 0, 0,   (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx))/2,          0, 0,  (Ls^2*(3*Fy + Ls*wy))/6,  1, 0;...
            0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  0, 1;...
            (wy*(Fx + Ls*wx)^4)/24, 0, 0, ((Fx + Ls*wx)^3*(Fy + Ls*wy))/6, 0, 0, (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx)^2)/4, 0, 0, (Ls^2*(3*Fy + Ls*wy)*(Fx + Ls*wx))/6, Fx + Ls*wx, 0, (Ls^3*(4*Fy + Ls*wy))/24, Ls, 1];
    end
    %% Solve for Consts Matrix
    
    Consts = pinv(Var_matrix)*bc_matrix;
    
end %end of case if index is empty or not





