function [sample,rf,rferr,inf,inferr,funcorr,funcorrerr,minfiberlength,maxfiberlength,positionall,freqall,repall,valuecorrall,eyesall,fibersall,etedall]=calculateparameters_withfiberid(file,globalallexDcut,globalallnum_pieces,globalalllength_pieces,thre1,thre2,thre3,interval,maxlength_bineyes,num_bineyes,limit_lengthfiber,fiber_id)

sample=[];
for i=1:length(file)
    
%The fuction calculatestruct.m calculate all the characteristic of the fibers and save them in the structure     
[globalallexDcut.(['exDcut' file{i}])]=calculatestruc(globalallexDcut.(['exDcut' file{i}]),globalallnum_pieces.(['num_pieces' file{i}]),thre1/(globalallexDcut.(['exDcut' file{i}])(1).unit_block),thre2/(globalallexDcut.(['exDcut' file{i}])(1).unit_block));

% I calculate the total replicated fraction, frequency of initiation and fork density for the time points
%"cellfun" apply the given function to each cell in cell array; the
%structure is coverted in cell array with the braket {}.
%The replicated fraction  for the experimental time point is calculated like the sum of all tha 1-blocks of
%all the fibers of a given experimantal time point divided by the sum of
%the lengths of the fibers.

lengthfiber(i)=sum(cellfun(@length,{globalallexDcut.(['exDcut' file{i}]).fiber}));
sumfiber(i)=sum(cellfun(@sum,{globalallexDcut.(['exDcut' file{i}]).fiber}));
exfrac(i)=sumfiber(i)/lengthfiber(i);
[globalallexDcut.(['exDcut' file{i}]).totrep_frac]=deal(exfrac(i));

for g=1:globalallnum_pieces.(['num_pieces' file{i}])
    globalallexDcut.(['exDcut' file{i}])(g).position=g;
    globalallexDcut.(['exDcut' file{i}])(g).file=file{i};
    globalallexDcut.(['exDcut' file{i}])(g).fiber_id=fiber_id{i}{g};
end

sample=[sample,globalallexDcut.(['exDcut' file{i}])(globalalllength_pieces.(['length_pieces' file{i}])>limit_lengthfiber)];
end


for i=1:length(sample)
% I calculate the total replicated fraction, frequency of initiation and fork density for the time points
%"cellfun" apply the given function to each cell in cell array; the
%structure is coverted in cell array with the braket {}.
%The replicated fraction is calculated like the sum of all tha 1-blocks of
%all the fibers of a given experimantal time point divided by the sum of
%the lengths of the fibers.
sample(i).lengthfiber=length(sample(i).fiber);
sample(i).sumfiber=sum(sample(i).fiber);

%I calculate the new origins on all the fibers of a time point.
%In the vectors containing the lengths of eyes for each fiber, "cellfun" recognize the eye lengths
%less than thre3 (converted in kb) and put a 1 at that positions and 0 at
%the others position. my is a cell array, with a cell for each fiber. 
% 'UniformOutput' in cell fun allow to combine outputs of different size
my=sample(i).length_eyes<=(thre3/(globalallexDcut.(['exDcut' file{1}])(1).unit_block));
sample(i).freq_init=sum(my)/((sample(i).lengthfiber-sample(i).sumfiber)*interval);

end


%I divide the fibers depending on replicated fraction and I calculate
%avereage replicated fraction, frequency of initiation and correlation
%function with errors

maxfiberlength=max([sample.lengthfiber]);
minfiberlength=min([sample.lengthfiber]);

centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);
reptot=[];
k=1;
delta=maxlength_bineyes/(num_bineyes);
for i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes
    select=[];
    position=[];
    select_freq=[];
    select_rep=[];
    valuecorrmatrix=zeros(1,2*maxfiberlength-1);
    for j=1:length(sample)
    temp=find(sample(j).fracfiber>i & sample(j).fracfiber<=i+(maxlength_bineyes/(num_bineyes)));
    select= [select,temp]; %Would select the time of the fiber; not usefull here, but in simulated fibers yes
    position= [position,j*ones(1,length(temp))];  %Select the position in sample
    sample(j).freq_init(isnan(sample(j).freq_init))=0;
    select_freq=[select_freq,sample(j).freq_init(temp)];
    sample(j).fracfiber(isnan(sample(j).fracfiber))=0;
    select_rep=[select_rep,sample(j).fracfiber(temp)];
    end
    for q=1:length(select)
    [valuecorrmatrix(q,maxfiberlength-sample(position(q)).lengthfiber+1:maxfiberlength+sample(position(q)).lengthfiber-1),lags]=xcorr(sample(position(q)).fiber',sample(position(q)).fiber','unbiased');
    end
        
    selectall{k}=select;
    positionall{k}=position;
    freqall{k}=select_freq;
    repall{k}=select_rep;
    valuecorrall{k}=valuecorrmatrix;
    valuecorrnozeros{k}=valuecorrmatrix(:,maxfiberlength-minfiberlength+1:maxfiberlength+minfiberlength-1);
    for z=1:length(positionall{k})
        eyesall{k}{z}=sample(positionall{k}(z)).length_eyes';
        etedall{k}{z}=sample(positionall{k}(z)).etedist';
        fibersall{k}{z}=sample(positionall{k}(z)).fiber;
    end
    k=k+1;
end



rf=cellfun(@mean,repall);
rferr=cellfun(@std,repall)./sqrt(cellfun(@length,repall));
inf=cellfun(@mean,freqall);
inferr=cellfun(@std,freqall)./sqrt(cellfun(@length,freqall));

for q=1:length(valuecorrall)
    for n=1:2*maxfiberlength-1
        if isempty(valuecorrall{q}(:,n))
            funcorr(q,n)=0;
            funcorrerr(q,n)=0;
        else
            funcorr(q,n)=mean(valuecorrall{q}(:,n));
            funcorrerr(q,n)=std(valuecorrall{q}(:,n))/sqrt(length(valuecorrall{q}(:,n)));
        end
    end
end


end

