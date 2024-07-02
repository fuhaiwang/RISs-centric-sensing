close all; clear
%% data load RIS2
load('RIS2_data.mat');
gain_final_all((gain_final_all>-10)==1) = []; 
A_all((gain_final_all>-10)==1) = []; 
%% initialization
gain_final_all_power = db2mag(gain_final_all+20);
[WT2,Ei,M,num_x,num_y] = RIS2_1819(A_all);
y1 = sqrt(gain_final_all_power);
E = Ei;     
algorithm = 'RWF';

b1 = WT2 * E;
b1_abs = abs(b1);
A = WT2;

%% data load
load('RIS1_data.mat');
gain_final_all((gain_final_all>-10)==1) = []; A_all((gain_final_all>-10)==1) = []; 
%% initialization
gain_final_all_power = sqrt(db2mag(gain_final_all+20));
[WT1,Ei,M] = RIS1_1819(A_all);
y2 = gain_final_all_power;
E = Ei;
o_y2 = abs(y2);

b2 = WT1 * E;
b2_abs = abs(b2);

%% STACK %% phasesolver 
% Run the Phase retrieval Algorithm
opts = struct;                  % Create an empty struct to store options
opts.algorithm = algorithm;     % Use the Fienup method to solve the retrieval problem.  Try changing this to 'twf' for truncated Wirtinger flow. 
opts.initMethod = 'optimal';    % Use the optimal spectral initializer method to generate an initial starting point for the solver  
opts.tol = 1e-2;                % The tolerance - make this smaller for more accurate solutions, or larger for faster runtimes
opts.verbose = 2;               % Print out lots of information as the solver runs (set this to 1 or 0 for less output)

fprintf('Running %s algorithm\n',opts.algorithm);
[X_, outs_stack_] = solvePhaseRetrieval([WT2;WT1], [], [b1_abs;b2_abs], M, opts);
pp_ = abs(X_)/max(abs(X_));
map_stack_simulated = rot90(reshape(pp_, [num_x,num_y]));
fprintf('Signal recovery required %d iterations (%f secs)\n',outs_stack_.iterationCount, outs_stack_.solveTimes(end));

fprintf('Running %s algorithm\n',opts.algorithm);
[X, outs_stack] = solvePhaseRetrieval([WT2;WT1], [], [y1;y2], M, opts);
pp_x = abs(X)/max(abs(X));
map_stack_measured = rot90(reshape(pp_x, [num_x,num_y]));
fprintf('Signal recovery required %d iterations (%f secs)\n',outs_stack.iterationCount, outs_stack.solveTimes(end));

%%
x = -5.2:2:4.8; 
y =    0:2:10; 
figure
subplot(1,2,1);
imagesc(x,y,flipud(map_stack_simulated));
colorbar
set(gca, 'YDir', 'normal');
set(gca, 'FontSize', 18);
hold on;
x_star = -1.2;
y_star = 6; 
star_size = 15;
plot(x_star, y_star, 'p', 'MarkerSize', star_size, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

subplot(1,2,2);
imagesc(x,y,flipud(map_stack_measured));
colorbar
set(gca, 'YDir', 'normal');

set(gca, 'FontSize', 18);
set(gcf,'unit','centimeters','position',[20 15 28 9.8])

hold on;
plot(x_star, y_star, 'p', 'MarkerSize', star_size, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

save_path = '.\result_stack\';
ax = gcf;
str=strcat(save_path,'2 RISs simulation and measured data stacked','.pdf');
exportgraphics(ax, str,'ContentType','vector');




