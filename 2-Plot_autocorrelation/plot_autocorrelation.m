
%___________________________PLOT AUTOCORRELATION_____________________________
%
% With this program I calculate the autocorrelation of each fiber and I set the
% number and size of bins in order to obtain the average profile of
% autocorrelation for each bin.
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_path='output_demo';
legend_graph='control';

%General variables
unit=1000; %Define the number of base pair (bp) for each block of the genome
v=0.5; %speed in kb/min

%Gaps smaller then thre1(here in bp)are combined in the analysis of fibers
thre1=1000; 

%Eyes smaller then thre2(here in bp)are not considered in the anamysis of
%fiber
thre2=1000;

%Eyes smaller than thre3(here in bp) are considered as new origins.
thre3=3000;
interval=thre3/(2*v*unit); %Detectable initiation events can occur in this interval (in min) for speed in kb/min
save([sample_path '/interval_' sample_path '.mat'],'interval');

% Fibers smaller of the threshold are removed from the analysis to avoid
% biases due to fiber size
limit_lengthfiber=80; 

%Number of bins by which the progam divide the fibers according to
%the replicated fraction
num_bineyes=7;
save([sample_path '/num_bineyes_' sample_path '.mat'],'num_bineyes');

%Maximum replicated fraction to use for the analysis
maxlength_bineyes=0.75;

%I load the data to analyse, that were previously saved by the function
%'storeexperfiber'
load(['../1-Data_extraction/' sample_path '/globalallexDcut.mat']);
load(['../1-Data_extraction/' sample_path '/globalallnum_pieces.mat']);
load(['../1-Data_extraction/' sample_path '/globalalllength_pieces.mat']);
load(['../1-Data_extraction/' sample_path '/file.mat']);
load(['../1-Data_extraction/' sample_path '/fiber_id.mat']);

centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);
i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Analysis untreated sample %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sample,rf,rferr,inf,inferr,funcorr,funcorrerr,minfiberlength,maxfiberlength,positionall,freqall,repall,valuecorrall,eyesall,fibersall,etedall]=calculateparameters_withfiberid(file,globalallexDcut,globalallnum_pieces,globalalllength_pieces,thre1,thre2,thre3,interval,maxlength_bineyes,num_bineyes,limit_lengthfiber,fiber_id);
save([sample_path '/sample_' sample_path '.mat'],'sample');
save([sample_path '/rf_' sample_path '.mat'],'rf');
save([sample_path '/rferr_' sample_path '.mat'],'rferr');
save([sample_path '/inf_' sample_path '.mat'],'inf');
save([sample_path '/inferr_' sample_path '.mat'],'inferr');
save([sample_path '/funcorr_' sample_path '.mat'],'funcorr');
save([sample_path '/funcorrerr_' sample_path '.mat'],'funcorrerr');
save([sample_path '/minfiberlength_' sample_path '.mat'],'minfiberlength');
save([sample_path '/maxfiberlength_' sample_path '.mat'],'maxfiberlength');
save([sample_path '/positionall_' sample_path '.mat'],'positionall');
save([sample_path '/freqall_' sample_path '.mat'],'freqall');
save([sample_path '/repall_' sample_path '.mat'],'repall');
save([sample_path '/valuecorrall_' sample_path '.mat'],'valuecorrall');
save([sample_path '/eyesall_' sample_path '.mat'],'eyesall');
save([sample_path '/fibersall_' sample_path '.mat'],'fibersall');
save([sample_path '/etedall_' sample_path '.mat'],'etedall');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
figure;
axis([0 0.75 0 Inf]) %  [left bottom width height]
hold on
ylabel('I(f) (1/(kb*min))','fontsize',12,'FontName','Arial');
xlabel('f','fontsize',12,'FontName','Arial');
errorbar(rf,inf,inferr,'LineStyle','none','Marker','o','MarkerSize',2.5,'linewidth',1);
legend(legend_graph);
im_paper1([sample_path '/' sample_path '_I'],4,3.3)


figure;
hold on
axis([0 80 0 0.75]) %  [left bottom width height]
ylabel('C(r,f)','fontsize',12,'FontName','Arial');
xlabel('r (kb)','fontsize',12,'FontName','Arial');
for i=1:length(funcorr(:,1))
x=-maxfiberlength+1:maxfiberlength-1;
H=errorbar(x,funcorr(i,:),funcorrerr(i,:));
errorbar_tick(H,1,'units')
end
legend([H],legend_graph);
im_paper1([sample_path '/' sample_path '_C'],4,3.3)

close all