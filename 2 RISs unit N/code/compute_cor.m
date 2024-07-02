function [risPosition,risRowVector,risColumnVector,RxPosition] = compute_cor(K, delta_theta)

    SOI_center_cor = [0;0;0];
    l = 25;
    delta_theta = 11.25;
    K = 16; 
    
    p_index = 0:K-1;
    p_theta = p_index*delta_theta; 
    
    
    risPosition = zeros(3, K);
    risRowVector = zeros(3, K);
    risColumnVector = zeros(3, K);
    RxPosition  = zeros(3, K);
    
    for i = 1: K
        risPosition_x_tmp = l * cosd(180+p_theta(i));
        risPosition_y_tmp = l * sind(180+p_theta(i));
        risPosition(:,i) = SOI_center_cor + [risPosition_x_tmp;risPosition_y_tmp;0];
    
        risRowVector(:,i) = [sind(180+p_theta(i)); -cosd(180+p_theta(i));0];
        risColumnVector(:,i) = [0;0;1];
    
    
        RxPosition_x_tmp = (l-1) * cosd(180+p_theta(i));
        RxPosition_y_tmp = (l-1) * sind(180+p_theta(i));
        RxPosition(:,i) = SOI_center_cor + [RxPosition_x_tmp;RxPosition_y_tmp;0];
    end



