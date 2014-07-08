function [F_matrix] = aif_exper_forces(Tags_bound_cond, Const_mat_combined , ws)

%first determine which tags forces will be returned for.
%F_matrix will also be a cellmat, holding forces Fx, Fy for each tag in
%each respective cell:
% F_matrix = {F_T1},{F_T2},{F_T3},{F_T4}



F_matrix = cellmat;


tag_data_length = length(Const_mat_combined);

for i = 1:tag_data_length
    
    
    if length(Const_mat_combined{i}) ~= 1
        
        Fx = [];
        Fy = [];
        Ls = Tags_bound_cond{1,2}(i); %Length of particular tag coil, this value can be used to better approximate force for terms in variable matrix that are coupled with Ls
        Tag_of_interest = i;
        for k = 1:size(Tags_bound_cond,1)
            
            %Find matrix for Boundary Conditions from tags
            bc_matrix = [0;... dy/dx(0), angle of center tag defined to be 0
                Tags_bound_cond{k,1}(1+Tag_of_interest);...  dy/dx(L), angle of Tag w.r.t. center tag (in center tag frame)
                0;... y(o), position of center tag, defined to be 0
                Tags_bound_cond{k,4}(Tag_of_interest)]; %y(L), position of Tag of interest
            
            
            %Find best approximation to variable matrix
            Var_matrix = bc_matrix * pinv(Const_mat_combined{i}); %Approximation to variable matrix
            
            
            %% Find Contribution from weight of spring, (not necessary for horizontal case, ws = 0
            coil_weight_den = ws*1e-3; %since ws is mN/mm and L is mm, then ws*L*1e-3 is Newton weight of spring, we will not multiply by length here as its accounted for in Bern, Eul. eqns
            %Now it is assumed that the weight is pulling down in the frame.
            
            wx =  -coil_weight_den*sin(-Tags_bound_cond{k,1}(1));
            wy = -coil_weight_den*cos(-Tags_bound_cond{k,1}(1));
            
            
            %% Find Fx and Fy from components of Var_Matrix
            %Use components of Variable matrix to find best
            %approximations to force
            
            Fy_2_13 = (1/3)*((Var_matrix(2,13)*6)/Ls^2 - Ls*wy);
            
            Fy_4_13 = (1/4)*((Var_matrix(4,13)*24)/(Ls^3) - Ls*wy);
            
%             Fy(end+1,1) = mean([Fy_2_13, Fy_4_13]);
            Fy(end+1,1) = Fy_2_13
            
            Fx_4_11 = Var_matrix(4,11)-Ls*wx; %possibly bad
            
            Fx_2_10 = (Var_matrix(2,10)*2)/(Ls*(2*Fy(end)+Ls*wy)) - Ls*wx;
            
            Fx_4_10 = (Var_matrix(4,10)*6)/(Ls^2*(3*Fy(end)+ Ls*wy)) - Ls*wx;
            
            
%             Fx(end+1,1) = mean([Fx_4_11, Fx_2_10, Fx_4_10]);
            Fx(end+1,1) = Fx_4_11;
            
            %         Var_matrix = [ 0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  1, 0;...
            %                        0, 0, 0,           (wy*(Fx + Ls*wx)^3)/6, 0, 0,      ((Fx + Ls*wx)^2*(Fy + Ls*wy))/2, 0, 0,   (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx))/2,          0, 0,  (Ls^2*(3*Fy + Ls*wy))/6,  1, 0;...
            %                        0, 0, 0,                               0, 0, 0,                                    0, 0, 0,                                    0,          0, 0,                        0,  0, 1;...
            %             (wy*(Fx + Ls*wx)^4)/24, 0, 0, ((Fx + Ls*wx)^3*(Fy + Ls*wy))/6, 0, 0, (Ls*(2*Fy + Ls*wy)*(Fx + Ls*wx)^2)/4, 0, 0, (Ls^2*(3*Fy + Ls*wy)*(Fx + Ls*wx))/6, Fx + Ls*wx, 0, (Ls^3*(4*Fy + Ls*wy))/24, Ls, 1];
            
        end
        
        F_matrix{i} = [Fx,Fy]; %where Fx and Fy or column vectors
        
    end
    
    
end

















