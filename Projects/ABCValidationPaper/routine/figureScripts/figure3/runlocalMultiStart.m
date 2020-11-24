R = ABCAddPaths('Rat_NPD','rat_STN_GPe');
%% Set Routine Pars
R.projectn = 'Rat_NPD'; % Project Name
R.out.tag = 'STN_GPe_ModComp'; % Task tag
R = simannealsetup_NPD_STN_GPe(R);

% Create Job
parallel.defaultClusterProfile('local')
c = parcluster();
j = createJob(c);
for i = 1:6
    createTask(j, @Figure3_i_rat_STN_GPe_MultiStartTest,1,{R});
end
% createTask(j, @sum, 1, {sum(randn(1,1e8))});
% createTask(j, @sum, 1, {sum(randn(1,1e8))});
% createTask(j, @sum, 1, {sum(randn(1,1e8))});

submit(j);
wait(j);
results = fetchOutputs(j)
delete(j);

