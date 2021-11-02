%% %% %% %% %% %% %% %% %% %% MASTER SCRIPT %% %% %% %% %% %% %% %% %% %% % %%
%%%                                                                        %%%
%%%    INFERENCE OF BRAIN NETWORKS WITH APPROXIMATE BAYESIAN COMPUTATION – %%%
%%%         FACE VALIDITY AND AN EXAMPLE APPLICATION IN PARKINSONISM       %%%
%%%                                                                        %%%
%%% Timothy O. West*, Luc Berthouze, Simon F. Farmer, Hayriye Cagnan,      %%%
%%% Vladimir Litvak                                                        %%%
%%%                                                                        %%%
%%% This script will reproduce the figures of the paper.                   %%%
%%% This software uses a number of external toolboxes, to whose authors we %%%
%%% are hugely grateful. The licenses for each are within the              %%%
%%% 'ABC_dependencies' folders of this repository.                         %%%
%%%                                                                        %%%
%%% The remainder of code is shared under a MIT License. For more details  %%%
%%% see the license file within this repository.                           %%%
%%%                                                                        %%%
%%% Timothy West, University of Oxford, 2020.                              %%%
%%% timothy.west@ndcn.ox.ac.uk                                             %%%
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%% %%

% This should point to your github repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox'; % your repo path
addpath(repopath)
addpath(genpath([repopath '\Projects\ABCValidationPaper\routine\figureScripts']))

% Setup the basic config structure 'R'
projname = 'ABCValidationPaper'; % This is the name of the current project
R = ABCAddPaths(repopath,projname);

%% Figure 2- 
% "Examining the convergence of ABC optimization upon summary statistics
% from recordings of the STN and GPe in Parkinsonian rats"
Figure2_i_NPD_ModelFitDemo(R);              % Model Fit Demo with STN/GPe
Figure2_ii_rat_STN_GPE_ModelComp(R);        % Model Comparison Demo with STN/GPe

%% Figure 3
%" Multi-start analysis to test face validity of the ABC-based estimation 
% of model parameters by demonstrating consistency of estimation and the 
% data specificity of parameter estimates."
Figure3_i_rat_STN_GPe_MultiStartTest(R);    % This will run the multistarts
Figure3_ii_analyseMultiStart(R);            % Multistart analysis            
Figure3_iii_plotMultiStart(R);              % Multistart plots
Figure3_iv_MultiStartStats(R);              % Multistart statistics

%% Figure 4
% "Testing face validity of the ABC model comparison approach to model
% identification." 
Figure4_i_rat_8node_Validation_ConfusionMat(R);
Figure4_ii_rat_8node_Validation_ConfusionMat_Evidence(R)

%% Figure 6
% "Scaling up the ABC model comparison framework – investigating models of
% the cortico-basal ganglia-thalamic network. "
Figure5_6_i_modelfitting(R)
Figure5_6_ii_modelcomparison

%% Appendix II
% "Examination of Integration Step-size"
FigureSI_testIntegrationStepSize

%% Appendix III
% "Examination of forward uncertainty of posterior model"
FigureSI_forwardUncertainty(R)