%% Setting up the COBRA toolbox and the GUROBI solver

rng(12) %this sets a seed value
% setting a seed value ensure that when generating random numbers 
% you are always generating the same random numbers
% this makes your results reproducible

% Initiate the gurobi solver 
cd /opt/gurobi810/linux64/matlab
gurobi_setup

% Initiate the Cobra toolbox
cd /home/helena/cobratoolbox/
initCobraToolbox

% Change the solver in to Cobra toolbox to gurobi
changeCobraSolver('gurobi','all')
%testAll %allows you to test the toolbox

% Change path to where the model and other matlab function
% files are saved 
cd /home/helena/PhD/Year4/Armida/

%% Importing the Arnold & Nikoloski 2014 leaf model
AraModel = readCbModel('Arabidopsis_Metabolism.xml'); 
% You can import any SBML model in this way

%% Extracting Model Features 
num2str(length(AraModel.rxns)) %no. of reactions
% it may be worth printing all of the reactions to get an 
% understanding of the model the kind of constraints that you can set
num2str(length(AraModel.mets)) %no. of metabolites

%% Setting experimental flux constraints 
% Setting and upper bound
AraModel = changeRxnBounds(AraModel,'Bio_opt',0.0005,'l');
% Setting a lower bound 
AraModel = changeRxnBounds(AraModel,'Bio_opt',0.001,'u');
% If you are setting an outflux for biomass, then you are parameterizing
% your model over the lifespan of the plant (up to the point where you 
% collected the biomass for measuring) 

% Other rxns that you might want to set constraints for:
% Im_NH4 (nitrogen uptake)
% RBC_h (CO2 assimilation)
% Bio_Nlim (biomass under nitrogen limited conditions) 

% You cannot combine metabolite constraint and biomass constraints in 
% the same model (because the time-scales are different)
% You can use total nitrogen in the leaf as a flux input constraint and
% combine that with the biomass constraints to see how it affects the
% model.

%% Doing a FBA Analysis
% checkObjective(AraModel); % checks what is the objective function

%Setting a new objective function (e.g. maximum CO2 assimilation)
AraModel = changeObjective (AraModel, 'RBC_h');
FBAsolution = optimizeCbModel(AraModel,'max');
FBAsolution.f % Maximum feasible value for the chosen objective function
FBAvalues = FBAsolution.v;
% You may want to save FBA results to a csv file to analyse your results
% FBAresults = [AraModel.rxns,printRxnFormula(AraModel, 'rxnAbbrList', AraModel.rxns),num2cell(FBAvalues)];
% writecell(FBAresults,'FBAResults.csv')

%% Doing a FVA Analysis
[minF, maxF] = fluxVariability(AraModel,95);
% Allows 95% of the optimal solution for the objective function
% You might want to check out some of the other setting in this function
% (for example you can set it to calculate with no loops)
FVAresults = [AraModel.rxns,num2cell(minF),num2cell(maxF)];
writecell(FVAresults,'FVAResults.csv')

%% Doing Flux Sampling 

iterations = 100; %no. of samples generates 
sDens = 10; %thinning

options.nStepsPerPoint = sDens; %sampling density 
options.nPointsReturned = iterations; %number of points returned 
options.toRound = 0; %whether or not the polytope is rounded
options.optPercentage = 95; %if 0 then no objective f'n is applied
[~, samples] =  sampleCbModel(AraModel, [], [], options); 
t_samples = transpose(samples);
csvwrite('SamplingResults.csv',t_samples)
