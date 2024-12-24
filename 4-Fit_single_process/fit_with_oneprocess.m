%_________________________FIT WITH ONE PROCESS___________________________
%
%With this program I do the fit with the genetic algorithm
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_path='output_demo';
num_of_fits=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(['../2-Plot_autocorrelation/' sample_path '/num_bineyes_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/interval_' sample_path '.mat']);

load(['../2-Plot_autocorrelation/' sample_path '/rf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/rferr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inferr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorrerr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/maxfiberlength_' sample_path '.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GENETIC ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% matlabpool open local 2
%Number of variables
nvars=3;

%Lower bounds
LB=[1e-10 1e-15 0];

%Upper bounds
UB=[10 1 5];

%Positions of integer variables
% IntCon=[1 2 4 8];

%Set options
%ga ignores or overwrites some options for mixed integer optimization problems
options = gaoptimset(@ga);

% Handle to the function that creates the initial population
% @gacreationuniform | @gacreationlinearfeasible
% options.CreationFcn=;

% Handle to the function that the algorithm uses to create crossover children
% @crossoverheuristic | {@crossoverscattered} | @crossoverintermediate | @crossoversinglepoint | @crossovertwopoint | @crossoverarithmetic
% options.CrossoverFcn={@crossoverintermediate, 1.5*ones(1,nvars)};

%The fraction of the population at the next generation, not including elite
%children, that is created by the crossover function (default=0.8)
options.CrossoverFraction=0.6;

%Positive integer specifying how many individuals in the current generation
%are guaranteed to survive to the next generation ({0.05*(default population size)} for mixed integer problems)
options.EliteCount=1;

%Scalar. If the fitness function attains the value of FitnessLimit, the
%algorithm halts. (default=-Inf)
% options.FitnessLimit=0.5;

%Hibrid faction
% hybridopts = optimset('Display','iter','MaxFunEvals',50);
% options.HybridFcn={@fmincon,hybridopts};

%Handle to the function that scales the values of the fitness function
% @fitscalingshiftlinear | @fitscalingprop | @fitscalingtop | default:@fitscalingrank
% options.FitnessScalingFcn=@fitscalingprop;

%Positive integer specifying the maximum number of iterations before the
%algorithm halts (default=100)
options.Generations=3000;

% Handle to a function that continues the optimization after ga terminates
% Function handle | @fminsearch | @patternsearch | @fminunc | @fmincon | {[]}
% options.HybridFcn=@fminsearch;

%Size of the population (default={min(max(10*nvars,40),100)} for mixed integer problems)
options.PopulationSize=[10*ones(1,10)];


%Direction of migration
% 'both' | defalut:'forward'
% options.MigrationDirection;

% Scalar between 0 and 1 specifying the fraction of individuals in each subpopulation that migrates to a different subpopulation
% default=0.2
options.MigrationFraction=0.1;

% Positive integer specifying the number of generations that take place between migrations of individuals between subpopulations
% default=20
options.MigrationInterval=50;

% Handle to the function that produces mutation children
% {@mutationuniform,rate] | @mutationadaptfeasible | {@mutationgaussian}
options.MutationFcn={@mutationadaptfeasible,0.3};

%Array of handles to functions that plot data computed by the algorithm
%@gaplotbestf | @gaplotbestindiv | @gaplotdistance | @gaplotexpectation | @gaplotgeneology
% | @gaplotmaxconstr | @gaplotrange | @gaplotselection | @gaplotscorediversity | @gaplotscores | @gaplotstopping 
options.PlotFcns={@gaplotbestf, @gaplotbestindiv, @gaplotscores, @gaplotdistance};

%Positive integer specifying the number of generations between consecutive
%calls to the plot functions (default=1)
% options.PlotInterval;

% Matrix or vector specifying the range of the individuals in the initial population 
options.PopInitRange=[LB;UB];


% Handle to the function that selects parents of crossover and mutation children
% @selectionremainder | @selectionuniform | {@selectionstochunif} |
% @selectionroulette | @selectiontournament 
options.SelectionFcn=@selectionroulette;

%Positive integer. The algorithm stops if there is no improvement in the objective function for StallGenLimit consecutive generations.
% default=50
options.StallGenLimit=3000;

% Positive scalar. The algorithm stops if there is no improvement in the objective function for StallTimeLimit seconds.
% {default=Inf}
% options.StallTimeLimit;

% Positive scalar. The algorithm stops after running for TimeLimit seconds.
% {default:Inf}
% options.TimeLimit;

% Positive scalar. TolCon is used to determine the feasibility with respect to nonlinear constraints.
% Positive scalar | {1e-6}
% options.TolCon=0.001;

% Positive scalar. The algorithm runs until the cumulative change in the fitness function value over StallGenLimit is less than TolFun.
% {default:1e-6}
%  options.TolFun=0.1;
 
% Compute fitness functions of a population in parallel.
% 'always' | {default='never'}
 options.UseParallel='always';


% String specifying whether the computation of the fitness function is
% vectorized (not possible in my case)
% 'on' | {default='off'} 
% options.Vectorized;




for i=[1:num_of_fits]
sprintf('Fit n. %i',i)
% Set the random number generetor for reproducibility
rng('shuffle');

tic,
clear garesult
[garesult.x,garesult.fval,garesult.exitflag,garesult.output,garesult.population,garesult.scores] = ga(@(var)functomin_oneprocess(rf,rferr,inf,inferr,funcorr,funcorrerr,maxfiberlength,var,interval),nvars,[],[],[],[],LB,UB,[],options);
toc

results(i,:)=[garesult.fval,garesult.x];

save([sample_path '/results_' sample_path '.mat'],'results');


end

