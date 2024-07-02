function [s_hat] = phase_retrival(abs_y, H, SOI_Num)
    PowofsigRecvTK_syn = abs_y;
    HofRISTK_Syn = H;
    M = SOI_Num;
    %{ 
    CoordinateDescent TWF Kaczmarz WirtFlow  RWF AmplitudeFlow TAF RAF    PhaseMax       
    PhaseLamp 
    PhaseLift SketchyCGM  GerchbergSaxton Fienup
    %} 

    opts = struct;                  % Create an empty struct to store options
    opts.algorithm = 'CoordinateDescent';      % Use the Fienup method to solve the retrieval problem.  Try changing this to 'twf' for truncated Wirtinger flow. 
    opts.initMethod = 'optimal';    % Use the optimal spectral initializer method to generate an initial starting point for the solver  
    opts.tol = 1e-2;                % The tolerance - make this smaller for more accurate solutions, or larger for faster runtimes
    opts.verbose = 2;               % Print out lots of information as the solver runs (set this to 1 or 0 for less output)
    
    % Run the Phase retrieval Algorithm
    fprintf('Running %s algorithm\n',opts.algorithm);
    [s_hat, outs] = solvePhaseRetrieval(HofRISTK_Syn, [], PowofsigRecvTK_syn, M, opts);
    fprintf('Signal recovery required %d iterations (%f secs)\n',outs.iterationCount, outs.solveTimes(end));


end