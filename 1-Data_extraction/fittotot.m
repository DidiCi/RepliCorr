
%___________________________DATA_EXTRACTION_____________________________
%
% With this program it is possible to extract the data from the excel files 
% obtained from the analysis of combed fibers with ImageJ. 
% The script calculates several global parameters for each experimental time point.
% The script allow to define the thresholds of fluorescence intensity, gaps size and eye
% length size. 
% The files globalallexDcut.mat, globalallnum_pieces.mat, globalalllength_pieces.mat
% will be used later in the analysis.
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Path to excel data, can contain mutiple path
path={'../Data_demo/timepoint1/' ...
      '../Data_demo/timepoint2/' ...
      '../Data_demo/timepoint3/' ...
      '../Data_demo/timepoint4/' ...
      '../Data_demo/timepoint5/'};
sample_path='output_demo';

%The function "intensities_extraction.m" extract the data from excel files
[intensities,fiber_id]=intensities_extraction(path);
save([sample_path '/intensities.mat'],'intensities')
%load([sample_path '/intensities.mat']);

%Threshold of fluorescence intensity, one for each data point
thre_int=[30*ones(1,3),80*ones(1,2)];

%Microscope conversion parameters
Conv_Zeiss100=0.16; %1pixel=0.16micrometer
Convmicro_kb=2;

%General variables
unit=1000; %Define the number of base pair (bp) for each block of the genome

%Variables for extraction
file={'timepoint_label1','timepoint_label2','timepoint_label3','timepoint_label4','timepoint_label5'}; %Name of the files I want to use
pag=1; %Pag from where to start to take the data; correspond to the page in the excel files with the informations on all the fibers

%Gaps smaller then thre1(here in bp)are combined in the analysis of fibers
thre1=1000; 

%Eyes smaller then thre2(here in bp)are not considered in the analysis of
%fiber
thre2=1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Intensities_treatment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I apply the threshold for intensities, gaps and eyes and I convert in Kb
[globalallexDcut,globalallnum_pieces,globalalllength_pieces]=intensities_treatment(intensities,file,unit,thre1,thre2,thre_int,Convmicro_kb,Conv_Zeiss100);
save([sample_path '/globalallexDcut.mat'],'globalallexDcut') 
save([sample_path '/globalallnum_pieces.mat'],'globalallnum_pieces')
save([sample_path '/globalalllength_pieces.mat'],'globalalllength_pieces') 
save([sample_path '/fiber_id.mat'],'fiber_id') 