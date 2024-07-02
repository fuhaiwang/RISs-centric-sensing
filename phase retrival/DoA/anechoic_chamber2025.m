% localization  
close all
clear
%%
load('matlab20.mat');
xx = gain_final_all>-10;
gain_final_all(xx==1) = []; 
A_all(xx==1) = []; 

gain_final_all_power = db2mag(gain_final_all+20);

degree = 0.5;
in_degree_1 = -8.5;  
out_degree = 13.15;  
N = 32; M = 180/degree+1;

Ei = zeros(180/degree+1,1);
Ei(round((90+in_degree_1)/degree)+1) = 1;  

light_speed = 3e8;
Fc = 5.8e9;
lambda = light_speed/Fc;

theta_i_M = -90:degree:90;
u_i = sin(theta_i_M./180*pi);

interUnits = 0.025; 
n = 1:N;
pT = -(n-1)*interUnits;
pTu_i = pT'*u_i;
steer_vector_in = exp(-1j*2*pi*pTu_i/lambda);

theta_s = out_degree/180*pi;  
u_s = sin(theta_s);
pTu_s = pT'*u_s;
steer_vector_out = exp(-1j*2*pi*pTu_s/lambda);
numofSamples = size(A_all,1);

index_Sel = [1];
temp = 1;

HofRIST = zeros(numofSamples,M);

for t=1:numofSamples
    phaseShiftMatrix_t = cell2mat(A_all(t));
    if sum(phaseShiftMatrix_t(1,:)) == 16
        index_Sel(temp) = t;
        temp = temp +1 ;
    end
    eta_RIS = ones(N,1).*0.9956 ;  
    phaseShiftMatrix = diag(exp(-1j*pi*flip(phaseShiftMatrix_t(10, 33:64))));
    phas_shiftvalue_ = diag(phaseShiftMatrix);
    eta_RIS(real(phas_shiftvalue_)==-1) = 0.9106;
    eta_RIS = diag(eta_RIS);
    HofRIST(t,:) = exp(1j*20/180*pi).*exp(-1j*2*pi*3/lambda).*transpose(steer_vector_out)*eta_RIS*phaseShiftMatrix*steer_vector_in.*cos(theta_i_M./180*pi);
end

WT = HofRIST;
rank(WT)

y = gain_final_all_power;

E = Ei;
Y = abs(WT*E);

%% phasesolver 
b = abs(WT * E);
A = WT;

opts = struct;                  % Create an empty struct to store options
opts.algorithm = 'RWF';      % Use the Fienup method to solve the retrieval problem.  Try changing this to 'twf' for truncated Wirtinger flow. 
opts.initMethod = 'optimal';    % Use the optimal spectral initializer method to generate an initial starting point for the solver  
opts.tol = 1e-2;                % The tolerance - make this smaller for more accurate solutions, or larger for faster runtimes
opts.verbose = 2;               % Print out lots of information as the solver runs (set this to 1 or 0 for less output)

% Run the Phase retrieval Algorithm
fprintf('Running %s algorithm\n',opts.algorithm);
[x, outs] = solvePhaseRetrieval(A, [], b, M, opts);
fprintf('Signal recovery required %d iterations (%f secs)\n',outs.iterationCount, outs.solveTimes(end));

figure
subplot(2,1,1);

plot(theta_i_M,abs(x)/max(abs(x)), '--','Color','#000000','LineWidth',1.5)
xlabel('degree(°) haha');
set(gcf,'unit','centimeters','position',[3 5 40 20])
vertical_line_x_1 = theta_i_M(round((90+in_degree_1)/degree)+1); % 竖线的位置
hold on;
% 3
o_y = abs(y);
[x_, outs] = solvePhaseRetrieval(A, [], o_y , M, opts);
ax = gca;
set(ax, 'FontSize', 16);
plot(theta_i_M,abs(x_)/max(abs(x_)),'Color','#2864F0','LineWidth',1.5)
xlabel('Degree($^{\circ}$ )', 'FontSize', 22,'interpreter','latex','FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 22,'interpreter','latex','FontWeight', 'bold');

set(gcf,'unit','centimeters','position',[3 5 34 16])
grid on

vertical_line_x_1 = theta_i_M(round((90+in_degree_1)/degree)+1); % 竖线的位置
hold on;
plot([vertical_line_x_1, vertical_line_x_1], [0, 1], 'r-','LineWidth',1.2);
hold off;
legend('Simulated data','Measured data','Real incident angle: -8.5$^{\circ}$','interpreter','latex','FontSize',18);

xlim([-90 90]);
ylim();

%% 

load('matlab25.mat');
xx = gain_final_all>-10;
gain_final_all(xx==1) = []; 
A_all(xx==1) = []; 
%%
gain_final_all_power = db2mag(gain_final_all+20);

degree = 0.5;
in_degree_1 = -13.25;  
out_degree = 8.5;  
N = 32; M = 180/degree+1;

Ei = zeros(180/degree+1,1);
Ei(round((90+in_degree_1)/degree)+1) = 1;  

light_speed = 3e8;
Fc = 5.8e9;
lambda = light_speed/Fc;

theta_i_M = -90:degree:90;
u_i = sin(theta_i_M./180*pi);

%%
interUnits = 0.025; 
n = 1:N;
pT = -(n-1)*interUnits;
pTu_i = pT'*u_i;
steer_vector_in = exp(-1j*2*pi*pTu_i/lambda);

theta_s = out_degree/180*pi; 
u_s = sin(theta_s);
pTu_s = pT'*u_s;
steer_vector_out = exp(-1j*2*pi*pTu_s/lambda);
numofSamples = size(A_all,1);

index_Sel = [1];
temp = 1;

HofRIST = zeros(numofSamples,M);

for t=1:numofSamples
    phaseShiftMatrix_t = cell2mat(A_all(t));
    if sum(phaseShiftMatrix_t(1,:)) == 16
        index_Sel(temp) = t;
        temp = temp +1 ;
    end
    eta_RIS = ones(N,1).*0.9956 ;  
    phaseShiftMatrix = diag(exp(-1j*pi*flip(phaseShiftMatrix_t(10, 1:32))));
    phas_shiftvalue_ = diag(phaseShiftMatrix);
    eta_RIS(real(phas_shiftvalue_)==-1) = 0.9106;
    eta_RIS = diag(eta_RIS);
    HofRIST(t,:) = exp(1j*20/180*pi).*exp(-1j*2*pi*3/lambda).*transpose(steer_vector_out)*eta_RIS*phaseShiftMatrix*steer_vector_in;%.*cos(theta_i_M./180*pi);
end

WT = HofRIST;      
rank(WT)
y = gain_final_all_power;
E = Ei;
Y = abs(WT*E);

%% phasesolver 
b = abs(WT * E);
A = WT;
opts = struct;                  % Create an empty struct to store options
opts.algorithm = 'RWF';      % Use the Fienup method to solve the retrieval problem.  Try changing this to 'twf' for truncated Wirtinger flow. 
opts.initMethod = 'optimal';    % Use the optimal spectral initializer method to generate an initial starting point for the solver  
opts.tol = 1e-2;                % The tolerance - make this smaller for more accurate solutions, or larger for faster runtimes
opts.verbose = 2;               % Print out lots of information as the solver runs (set this to 1 or 0 for less output)

% Run the Phase retrieval Algorithm
fprintf('Running %s algorithm\n',opts.algorithm);
[x, outs] = solvePhaseRetrieval(A, [], b, M, opts);
fprintf('Signal recovery required %d iterations (%f secs)\n',outs.iterationCount, outs.solveTimes(end));

subplot(2,1,2);
plot(theta_i_M,abs(x)/max(abs(x)), '--','Color','#000000','LineWidth',1.5)

set(gcf,'unit','centimeters','position',[3 5 40 20])
vertical_line_x_1 = theta_i_M(round((90+in_degree_1)/degree)+1); % 竖线的位置
hold on;

o_y = abs(y);
[x_, outs] = solvePhaseRetrieval(A, [], o_y , M, opts);
ax = gca;
set(ax, 'FontSize', 16);
plot(theta_i_M,abs(x_)/max(abs(x_)),'Color','#2864F0','LineWidth',1.5)
xlabel('Degree($^{\circ}$ )', 'FontSize', 22,'interpreter','latex','FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 22,'interpreter','latex','FontWeight', 'bold');

set(gcf,'unit','centimeters','position',[3 5 25 20])
grid on

vertical_line_x_1 = theta_i_M(round((90+in_degree_1)/degree)+1); % 竖线的位置
hold on;
plot([vertical_line_x_1, vertical_line_x_1], [0, 1], 'r-','LineWidth',1.2);
hold off;
legend('Simulated data','Measured data','Real incident angle: -13.25$^{\circ}$','interpreter','latex','FontSize',18);

xlim([-90 90]);
ylim();

save_path = '.\result\';
ax = gcf;
str=strcat(save_path,'RIS Reconstructed simulation angle RIS12','.pdf');
exportgraphics(ax, str,'ContentType','vector');













