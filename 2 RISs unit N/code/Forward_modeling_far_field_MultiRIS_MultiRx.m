close all % N
clear
%% RIS
light_speed = 3e8; 
Fc = 30e9;
Lambda_c   = light_speed/Fc;
InterUnits = Lambda_c*1;
% SOI Tx 
SOI_Position= [-500;-500;0]./100;
SOI_XVector = [1;0;0];
SOI_YVector = [0;1;0];
SOI_ZVector = [0;0;1];
SOI_XNum    = 10;
SOI_YNum    = 10;
SOI_ZNum    = 1;
SOI_Num     = SOI_XNum*SOI_YNum*SOI_ZNum;
SOI_Resolution = 100./100;
Ei = zeros(SOI_Num, 1, 'double');  
[letter_index, letter_coor] = Font_library_3D_linear("star", SOI_Position, SOI_XVector, SOI_YVector, SOI_ZVector, SOI_XNum, SOI_YNum, SOI_ZNum, SOI_Resolution); 

Ei(letter_index) = [1,1,1,1,0.5,0.5,0.2,0.2,1,1,1,1,0.3,0.8,1,0.8,0.8,0.8,0.65,0.65,0.65]; 
TxPosition = letter_coor; 
%_________________________
risRowColNum = 4;
svd_all = zeros(risRowColNum,SOI_Num);
ls_matrix = zeros(risRowColNum,SOI_ZNum,SOI_XNum,SOI_YNum);
E_asyn_matrix = zeros(risRowColNum,SOI_ZNum,SOI_XNum,SOI_YNum);
E_syn_matrix = zeros(risRowColNum,SOI_ZNum,SOI_XNum,SOI_YNum);
delt_RISNum = 10;

for risrowcolNum = 1:risRowColNum
%%
% RIS
    close all
    X =sprintf('RISs block is %d',risrowcolNum);  
    disp(X)
    RISBolckNum = 4; 
    RxNum = 1; 
    K = RISBolckNum;
    
    [risPosition_all,risRowVector_all,risColumnVector_all,RxPosition_all] = compute_cor(16, 180/16); 
    
    risPosition = zeros(3,RISBolckNum);
    risRowVector = zeros(3,RISBolckNum);
    risColumnVector = zeros(3,RISBolckNum);
    RxPosition = zeros(3,RISBolckNum);

    for j = 1:RISBolckNum % 
        i_ = 1+(j-1)*16/RISBolckNum;
        risPosition(:,j) = risPosition_all(:,i_); 
        risRowVector(:,j) = risRowVector_all(:,i_); 
        risColumnVector(:,j) = risColumnVector_all(:,i_); 
        RxPosition(:,j) = RxPosition_all(:,i_); 
    end 

    risRowNum = 10 + risrowcolNum*delt_RISNum; 
 
    risColNum = 1;
    risNum = risRowNum*risColNum;
    for p = 1:RISBolckNum  
        xx = risRowVector(:,p)/norm(risRowVector(:,p));
        risPosition(:,p) = risPosition(:,p) - 1/2.*xx*risRowNum.*InterUnits;
    end
    % %% 3.Visualization
    Ei_hat = zeros(1,1);
    Visualize_MIMO_RIS(SOI_Position,SOI_XVector,SOI_YVector,SOI_ZVector,SOI_XNum,SOI_YNum,SOI_ZNum,SOI_Resolution,letter_coor,letter_index,Ei,risPosition(:,1:RISBolckNum),risRowVector(:,1:RISBolckNum),risColumnVector(:,1:RISBolckNum),risRowNum,risColNum,RxPosition(:,1:RISBolckNum),Ei_hat,risNum);
    %% modeling
    numofSamples = 100;

    HofRIST = zeros(numofSamples,SOI_Num);
    sigRecvT = zeros(numofSamples,1);
    
    HofRISTK = zeros(numofSamples*K,SOI_Num);
    HofRISTK_Syn = zeros(numofSamples,SOI_Num);
    sigRecvTK = zeros(numofSamples*K, 1);
    
    for k = 1:K
        X =sprintf('RIS %d is setting up and sampling...', k);
        disp(X)
        ris = ClassRISMISO_FarField_new(risRowNum, risColNum, Fc, risPosition(:,k), risRowVector(:,k), risColumnVector(:,k)); 
        units_position_ris = ris.getUnitPosition;
    
        [theta_i_M, phi_i_M] = ris.SOI_RIS_Incidence_angle(TxPosition);
        [theta_s, phi_s] = ris.SOI_RIS_Incidence_angle(RxPosition(:,k));
        for t = 1: numofSamples
            wRandom = randi([0,1],1,risNum);
            [sigRecv1, H_temp] = ris.takeReflection(RxNum, TxPosition, RxPosition(:,k), theta_i_M, phi_i_M, theta_s, phi_s, risPosition(:,k),units_position_ris, wRandom, Ei);
            
            HofRIST(t,:) = H_temp;
            sigRecvT(t,:) = sigRecv1;
        end
        HofRISTK((k-1)*numofSamples+1:k*numofSamples,:) = HofRIST;   % [RxNum*numofSamples, SOI_Num]
        sigRecvTK((k-1)*numofSamples+1:k*numofSamples,:) = sigRecvT; % [RxNum*numofSamples, 1]
        HofRISTK_Syn = HofRISTK_Syn + HofRIST; 
    end
    %% [U, S, V] = svd(A);
    [U_1, S_1, V_1] = svd(HofRISTK);
    [U_2, S_2, V_2] = svd(HofRISTK_Syn);
    svd_all(risrowcolNum,:) = diag(S_1);

    SNR_all = [30];
    SN = mean(Ei.^2);  
    SNR = SNR_all(1);
    sigma_Nosie = SN / 10^(SNR(1)/10);
    W = complex(randn(numofSamples*K,1), randn(numofSamples*K,1)); 
    
    sigRecvTK_asyn = HofRISTK*Ei;% + sqrt(sigma_Nosie/2)*W;

    sigRecvTK_syn = HofRISTK_Syn*Ei; 
    %% disp
    rank_HofRISTK = fprintf('M is %d,E is %d, and the rank of the measurement operator AW is %d',SOI_Num,size(letter_index,2),rank(HofRISTK));
    disp(rank_HofRISTK);
    rank_HofRISTKN(risrowcolNum) = rank(HofRISTK,1e-11);
    %% solving 
    WT = HofRISTK;
    y = sigRecvTK_asyn;
    E = Ei;
    Mx = SOI_XNum;
    My = SOI_YNum;
    Mz = SOI_ZNum;
    E_x_ls = pinv(WT)*y;
    E_x_ls_abs = abs(E_x_ls); 
    %% Criterion RE（dB）AND SNR 
    RE(risrowcolNum) = mean((E-E_x_ls_abs).^2)/mean(E.^2);
    sigma2 = 1;
    SNR(risrowcolNum) = 10*log10(mean(E.^2)/sigma2);
    %% plot and save
    save_path = '..\result\';
    for layer = 1:Mz
        E_matrix = reshape(E(1+Mx*My*(layer-1):Mx*My*layer), [Mx, My])';
        E_x_ls_abs_matrix = reshape(E_x_ls_abs(1+Mx*My*(layer-1):Mx*My*layer), [Mx, My])';
        
        ls_matrix(risrowcolNum,layer,:,:)= flipud(E_x_ls_abs_matrix);
        ax = gcf;
        str=strcat(save_path,'N=_',num2str(risrowcolNum),'_',num2str(risRowNum),'Far field layer',num2str(layer),'.pdf');
        SSIM(risrowcolNum,layer) = ssim(E_x_ls_abs_matrix,E_matrix);
    end
end
%%
figure
colors = {'b', 'g', 'm', 'c', 'r', 'y', 'k', [0.5 0.5 0.5], [0.2 0.8 0.2], [0.8 0.4 0.6]};
linestyles = {'-.', '-.', '-.', '-.', '-', '-', '--', ':', '-.', '-', '--'};  
hLegend = [];
for i =1:size(svd_all,1)
    color = colors{i}; 
    linestyle = linestyles{i}; 
    hLine = plot(svd_all(i,:),'linestyle',linestyle,'LineWidth', 2, 'Color', color);
    hold on 
    hLegend = [hLegend, hLine];
end
index_ = 10+(1:4)*delt_RISNum;
legend(hLegend,['N=',num2str(index_(1))],['N=',num2str(index_(2))],['N=',num2str(index_(3))],['N=',num2str(index_(4))], 'FontSize', 16,'interpreter','latex');
ylabel('Singular value','FontSize',12,'interpreter','latex');
set(gca, 'FontSize', 18);

set(gca, 'FontSize', 18);
grid on
ax = gcf;
str=strcat(save_path,'N_Singular_value_distribution','.pdf');
exportgraphics(ax, str,'ContentType','vector');

%%
index = index_;
figure
yyaxis left;
plot(index,SSIM','--*','LineWidth', 2);
xlabel('Scenario index', 'FontSize', 12,'interpreter','latex');
ylabel('SSIM value', 'FontSize', 12,'interpreter','latex');
xlim([index(1)-0.2, index(end)+0.2]);  
ylim([0.8, 1.005]); 
set(gca, 'FontSize', 18);

grid on
yyaxis right;
plot(index, RE,'--o','LineWidth', 2);
legend('SSIM','Relative error', 'FontSize', 16,'interpreter','latex','Location','northeast');
xlabel('$N_k$','FontSize',12,'interpreter','latex');
ylabel('Relative error','FontSize',12,'interpreter','latex');
xlim([index(1)-0.2, index(end)+0.2]); 
ylim([-.01, 0.12]); 

set(gca, 'FontSize', 18);
grid on
ax = gcf;
str=strcat(save_path,'N_SSIM','.pdf');
exportgraphics(ax, str,'ContentType','vector');
%%
figure
for i=1:risRowColNum
    for j=1:SOI_ZNum
        subplot(1,risRowColNum,i);
        set(gcf,'unit','normalized','position', [0.05,0.05,0.5,0.18]);
        imagesc(squeeze(ls_matrix(i,j,:,:)));  
        colorbar
    end
    set(gca, 'FontSize', 15);
end

ax = gcf;
str=strcat(save_path,'N_ls_matrix','.pdf');
exportgraphics(ax, str,'ContentType','vector');












