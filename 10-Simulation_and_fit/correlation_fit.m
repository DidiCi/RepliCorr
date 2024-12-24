%I take the data from sample, I calculate the freq.of initiation and the
%correlation function in the indicated interval
%Then I fit the data with a genetic algorithm
%The speed is in kb/sec, the frquency of initiation in 1/(kb*sec)

%no automatic integer setting
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%Options of genetic algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of variables
nvars=3;

%Order of variables v,I0,alpha
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
options.FitnessLimit=0.5;

%Hibrid faction
% hybridopts = optimset('Display','iter','MaxFunEvals',50);
% options.HybridFcn={@fmincon,hybridopts};

%Handle to the function that scales the values of the fitness function
% @fitscalingshiftlinear | @fitscalingprop | @fitscalingtop | default:@fitscalingrank
% options.FitnessScalingFcn=@fitscalingprop;

%Positive integer specifying the maximum number of iterations before the
%algorithm halts (default=100)
options.Generations=1500;

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
options.StallGenLimit=1500;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load('sample.mat');

vel=sample(1).vel;
interval=sample(1).interval;
Fiber_number=sample(1).Fiber_number;
Fiber_length=sample(1).Fiber_length;
Init_prob=sample(1).Init_prob;
unit=sample(1).unit;
thre=sample(1).thre;

%%%%%%%Calculation of frequency of initiation and correlation%%%%%%%%%%%%

%I divide the fibers depending on replicated fraction and I calculate
%avereage replicated fraction, frequency of initiation and correlation
%function with errors
maxlength_bineyes=1;
num_bineyes=10;
centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);
reptot=[];
k=1;
delta=maxlength_bineyes/(num_bineyes)
for i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes
    select=[];
    position=[];
    select_freq=[];
    select_rep=[];
    valuecorrmatrix=[];
    realselect_freq=[];
    for j=1:Fiber_number
    temp=find(sample(j).rep_frac>i & sample(j).rep_frac<=i+(maxlength_bineyes/(num_bineyes)));
    select= [select,temp];
    position= [position,j*ones(1,length(temp))];
    sample(j).freq_init(isnan(sample(j).freq_init))=0;
    select_freq=[select_freq,sample(j).freq_init(temp)];
    sample(j).realfreq_init(isnan(sample(j).realfreq_init))=0;
    realselect_freq=[realselect_freq,sample(j).realfreq_init(temp)];
    sample(j).rep_frac(isnan(sample(j).rep_frac))=0;
    select_rep=[select_rep,sample(j).rep_frac(temp)];
    end
    for q=1:length(select)
    [valuecorrmatrix(q,:),lags]=xcorr(sample(position(q)).fiber(select(q),:),sample(position(q)).fiber(select(q),:),'unbiased');
    end
    selectall{k}=select;
    positionall{k}=position;
    freqall{k}=select_freq;
    realfreqall{k}=realselect_freq;
    repall{k}=select_rep;
    valuecorrall{k}=valuecorrmatrix;
    k=k+1;
end

rf=cellfun(@mean,repall);
rferr=cellfun(@std,repall)./sqrt(cellfun(@length,repall));
inf=cellfun(@mean,freqall);
inferr=cellfun(@std,freqall)./sqrt(cellfun(@length,freqall));
realinf=cellfun(@mean,realfreqall);
realinferr=cellfun(@std,realfreqall)./sqrt(cellfun(@length,realfreqall));
figure;
hold on
ylabel('I(t)');
xlabel('f(t)');
errorbar(rf,inf,inferr);

for q=1:length(valuecorrall)
   funcorr(q,:)=mean(valuecorrall{q});
   funcorrerr(q,:)=std(valuecorrall{q})/sqrt(length(valuecorrall{q}(:,1)));
end

figure;
hold on
x=-Fiber_length+1:Fiber_length-1;
ylabel('Corr');
xlabel('x');
for i=1:length(funcorr(:,1))
    plot(x,funcorr(i,:));
% errorbar(x,funcorr(i,:),funcorrerr(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('results.mat');
%%%%%%%%%%%%%%%%%%%%%%%%% Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=[1:100]
% % Set the random number generetor for reproducibility
rng('shuffle');

tic,
clear garesult
[garesult.x,garesult.fval,garesult.exitflag,garesult.output,garesult.population,garesult.scores] = ga(@(var2)functomin(rf,rferr,inf,inferr,funcorr,funcorrerr,Fiber_length,var2,interval),nvars,[],[],[],[],LB,UB,[],options);
toc

results(i,:)=[garesult.fval,garesult.x]
save('results.mat','results');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

