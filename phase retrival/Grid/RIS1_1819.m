function [WT,Ei,M] = RIS1_1819(A_all) 
    N = 32; 

    axis_x = -5.2:2:4.8;  
    ayis_y = 0:2:10;


    M = length(axis_x)*length(ayis_y);
    [xx, yy]= meshgrid(axis_x, ayis_y); 
    

    RIS1_cor = [1.405,0.02]; 
    n1 = [-0.11,0.6];
    resolution = 1;
    yyy = yy.*resolution-RIS1_cor(2);
    xxx = xx.*resolution-RIS1_cor(1);
    
    r_i = sqrt((xxx - RIS1_cor(1)).^2 + (yyy - RIS1_cor(2)).^2) ;
    r_i_vec = diag(reshape(r_i',[M,1]));

    theta_i_M_ = (atan2(yyy, xxx) - atan2(n1(2), n1(1)))*180/pi;
    theta_i_M = reshape(theta_i_M_', [1,length(axis_x)*length(ayis_y)]);
    
    y_ = yy.*resolution;
    x_ = xx.*resolution;

    y_index_matrix = min(abs(y_-6))==abs(y_-6);
    index_y = find(y_index_matrix(:,1), 1);
    x_index_matrix = min(abs(x_'+1.2))==abs(x_'+1.2);
    index_x = find(x_index_matrix(:,1), 1);

    in_degree_index = (index_y-1)*length(axis_x)+index_x;  
    out_degree = -8.5; 
    
    Ei = zeros(length(axis_x)*length(ayis_y),1);
    Ei(in_degree_index) = 1; 
    
    light_speed = 3e8;Fc = 5.8e9; 
    lambda = light_speed/Fc;
    u_i = sin(theta_i_M./180*pi);
    %% model
    interUnits = 0.025; 
    n = 1:N;
    pT = -(n-1)*interUnits;
    pTu_i = pT'*u_i;
    steer_vector_in = exp(-1j*2*pi*pTu_i/lambda);
    
    theta_s = out_degree/180*pi;  
    u_s = sin(theta_s);
    pTu_s = pT'.*u_s;
    steer_vector_out = exp(-1j*2*pi*pTu_s/lambda);
    numofSamples = size(A_all,1);
    
    index_Sel = [1];temp = 1;
    HofRIST = zeros(numofSamples, M);
    
    for t=1:numofSamples
        phaseShiftMatrix_t = cell2mat(A_all(t));
        phaseShiftMatrix_t = phaseShiftMatrix_t(10, 1:32);
        if sum(phaseShiftMatrix_t(1,:)) == 16 
            index_Sel(temp) = t;
            temp = temp +1 ;
        end
        eta_RIS = ones(N,1).*0.9956;    % 0=0.9956，pi=0.9106
        phaseShiftMatrix = diag(exp(-1j*pi*flip(phaseShiftMatrix_t(1,:))));
        eta_RIS(diag(real(phaseShiftMatrix))==-1) = 0.9106; 
        eta_RIS = diag(eta_RIS);
        HofRIST(t,:) = transpose(steer_vector_out)*eta_RIS*phaseShiftMatrix*steer_vector_in*exp(-1j*2*pi*r_i_vec/lambda);
    end

    WT = HofRIST;
end


















