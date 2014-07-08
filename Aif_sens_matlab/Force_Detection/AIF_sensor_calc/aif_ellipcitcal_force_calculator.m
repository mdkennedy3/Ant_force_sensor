function [Fmag, Fx, Fy] = aif_ellipcitcal_force_calculator(pos_prof, Const, Tag_of_interest)

% Consts should contain 6 elements, 
% pos_prof is the profile of positions to be used to find:
% Force Magnitude: Fmag
% Force x comp: Fx
% Force y comp: Fy

%% Setup:

%Set C's
C1 = Const(1);
C2 = Const(2);
C3 = Const(3);
C4 = Const(4);
C5 = Const(5);
C6 = Const(6);

F_mu = C3 + C6 - (C2^2)/(4*C1) - (C5^2)/(4*C4);

mu_x = (C2)/(2*C1);

mu_y = (C5)/(2*C4);

%% Find magnitudes for all positions provided with Const model:
Fmag = [];
for i = 1:size(pos_prof,2)
    if ~isempty(pos_prof{1+Tag_of_interest,i})
    xp =  pos_prof{1+Tag_of_interest,i}(1); %mm
    yp =  pos_prof{1+Tag_of_interest,i}(2); %mm
   Fmag(end+1,1) = C1 *(xp + mu_x)^2 + C4*(yp + mu_y)^2 + F_mu; 
    end
end

%% Find Force Vector for each position provided with Const model:
F_vec = [];
for i = 1:size(pos_prof,2)
    if ~isempty(pos_prof{1+Tag_of_interest,i})
    xp =  pos_prof{1+Tag_of_interest,i}(1); %mm
    yp =  pos_prof{1+Tag_of_interest,i}(2); %mm

F_vec(end+1,:) = (C1 *(xp + mu_x)^2 + C4*(yp + mu_y)^2 + F_mu)/(sqrt((2*C1*(xp + mu_x))^2 + (2*C4*(yp + mu_y))^2)).*[2*C1*(xp + mu_x) , 2*C4*(yp + mu_y) ]; %mN
    end
end

Fx = F_vec(:,1); %mN
Fy = F_vec(:,2); %mN