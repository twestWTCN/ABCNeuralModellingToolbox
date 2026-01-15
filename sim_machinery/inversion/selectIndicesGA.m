function optimalIndices = selectIndicesGA(R,parBank,pIndMap,pOrg)
% GA options - Adjust these according to your problem's specifics
options = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 150, ...
    'EliteCount', 2, 'CrossoverFraction', 0.8, ...
    'MutationFcn', @mutationuniform, 'Display', 'none'); %,'PlotFcn', {@gaplotdistance, @gaplotbestf});

% Binary representation for sample selection
% 1 means the sample is selected, 0 means it is not
numSamples = size(parBank,2);
IntCon = 1:numSamples; % Treat every variable as integer (binary in this case)
disp('Using GA to optimize sample draws for proposal!')
% Running GA
[optimalIndices, fval] = ga(@(indices) sampleFitness(indices,R,parBank,pIndMap,pOrg), numSamples, [], [], [], [], ...
    zeros(1, numSamples), ones(1, numSamples), [], IntCon, options);
