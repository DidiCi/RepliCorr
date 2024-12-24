%I simulate a nucleation and growth process on N blocks with variable
%velocity and I(t). 
%The initiation events are selected depending on the threshold thre. The
%I(t) for the fiber is then calculated as for the experiments.
%The unbiased autocorrelation and the I(t) are plotted and the set of
%simulated fibers saved in sample

clear all;
close all;

Fiber_length=150;
Fiber_number=100;

thre=3; %threshold for initiation
unit=1000; %unit of the block in bp
vel=1; %Num blocks elongated per round -> used in simulation;

%Reference for unit of time
interval=1; %%interval in min between two steps of simiulation
vel_permin=vel*unit/interval; %Speed express in bp/min
realinterval=thre*interval; %interval in min in which one event can be detected; 
coeff=0;
Init_prob=@(time) (0.03*interval^(coeff+1))*time.^coeff; %I0*t^alpha with alpha in 1/(block*interval^(alpha+1))==probability/interval

%Set flag=1 to generate set of fibers
flag=1;
if flag==1
% I replicate the fibers
for i=1:Fiber_number
   clear fiber rep_frac left_fork right_fork initiations eyes real_init rep_fracinit;
   step=1;
   Init=Init_prob(interval);
   fiber=zeros(1,Fiber_length);
   rep_frac(1)=sum(fiber(1,:))/Fiber_length; 
   temp=ceil(rand(1,round(Init*Fiber_length))*Fiber_length);
   initiations{1}=temp;
  %Find forks
  left_fork{1}=(findstr(fiber(1,:),[0 1])+1)';
  right_fork{1}=(findstr(fiber(1,:),[1 0]))';
 
  while sum(fiber(end,:))<Fiber_length
      Init=Init_prob(step*interval+interval);
      fiber(step+1,:)=fiber(step,:);
      %Elongation
      %I elongate the origins of the previous step
      for k=0:vel-1
      fiber(step+1,initiations{step}((initiations{step}-k)>=1)-k)=1;
      end
      for k=1:vel
      fiber(step+1,initiations{step}((initiations{step}+k)<=Fiber_length)+k)=1;
      end
      %Not problematic forks
      for k=1:vel
      fiber(step+1,left_fork{step}(2:end)-k)=1;
      fiber(step+1,right_fork{step}(1:end-1)+k)=1;         
      end
      %First left fork
      if ~isempty(left_fork{step})
      if  left_fork{step}(1)-vel<=0
         fiber(step+1,1:left_fork{step}(1))=1;
      else
          for k=1:vel
            fiber(step+1,left_fork{step}(1)-k)=1;     
          end
      end
      end
      %Last right fork
      if ~isempty(right_fork{step})
      if  right_fork{step}(end)+vel>=Fiber_length
         fiber(step+1,right_fork{step}(end):Fiber_length)=1;
      else
          for k=1:vel
            fiber(step+1,right_fork{step}(end)+k)=1;   
          end
      end  
      end
      %Initiation
      temp=find(rand(1,Fiber_length)<Init);
      new_init=[];
      for p=1:length(temp)
          if fiber(step+1,temp(p))==0
            new_init=[new_init,temp(p)];
          end
      end
%       temp=find(fiber(step,:)==0);
%       if length(temp)<=round(Init*length(temp))
%           new_init=temp;
%       else
%           new_init=temp(ceil(rand(1,round(Init*length(temp)))*length(temp)));
%       end
      initiations{step+1}=new_init; %The saved initiations are on unreplicated DNA
      %Find new forks  
      left_fork{step+1}=(findstr(fiber(step+1,:),[0 1])+1)';
      right_fork{step+1}=(findstr(fiber(step+1,:),[1 0]))';
      rep_frac(step+1)=sum(fiber(step+1,:))/Fiber_length;
      step=step+1;
  end
  %Calculate eyes length
  for j=1:length(fiber(:,1))
      if fiber(j,1)==0 && fiber(j,end)==0 , eyes{j}=right_fork{j}-left_fork{j}+1;
      elseif fiber(j,1)==1 && fiber(j,end)==0 && ~isempty(right_fork{j}(2:end)) , eyes{j}=right_fork{j}(2:end)-left_fork{j}+1;
      elseif fiber(j,1)==1 && fiber(j,end)==0 && isempty(right_fork{j}(2:end))  , eyes{j}=[];  
      elseif fiber(j,1)==0 && fiber(j,end)==1 && ~isempty(left_fork{j}(1:end-1))  , eyes{j}=right_fork{j}-left_fork{j}(1:end-1)+1;
      elseif fiber(j,1)==0 && fiber(j,end)==1 && isempty(left_fork{j}(1:end-1))  , eyes{j}=[]; 
      elseif fiber(j,1)==1 && fiber(j,end)==1 , eyes{j}=right_fork{j}(2:end)-left_fork{j}(1:end-1)+1;
      end
      real_init{j}=eyes{j}(eyes{j}<=thre);
  end
  
  sample(i).fiber=fiber;
  sample(i).rep_fracinit=rep_frac;
  sample(i).rep_frac=rep_frac;
  sample(i).left_fork=left_fork;
  sample(i).right_fork=right_fork;
  sample(i).eyes=eyes;
  sample(i).real_init=real_init;
  sample(i).initiations=initiations;
  sample(i).freq_init=cellfun(@length,sample(i).initiations)./((Fiber_length-sample(i).rep_frac*Fiber_length)*interval); 
  sample(i).realfreq_init=cellfun(@length,sample(i).real_init)./((Fiber_length-sample(i).rep_fracinit*Fiber_length)*realinterval);
end


sample(1).vel=vel;
sample(1).interval=interval;
sample(1).Fiber_number=Fiber_number;
sample(1).Fiber_length=Fiber_length;
sample(1).Init_prob=Init_prob;
sample(1).unit=unit;
sample(1).thre=thre;

save('sample.mat','sample');
end 

load('sample.mat');

%I plot the evolution of the fiber
figure
Num=1;
hold all
surface(1:Fiber_length,1:length(sample(Num).fiber(:,1)),sample(Num).fiber)
hold off

%I divide the fibers depending on replicated fraction and I calculate
%avereage replicated fraction, frequency of initiation and correlation
%function with errors
maxlength_bineyes=1;
num_bineyes=10;
centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);
reptot=[];
k=1;

for i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes
    select=[];
    position=[];
    select_freq=[];
    select_rep=[];
    valuecorrmatrix=[];
    for j=1:Fiber_number
    temp=find(sample(j).rep_frac>i & sample(j).rep_frac<=i+(maxlength_bineyes/(num_bineyes)));
    select= [select,temp];
    position= [position,j*ones(1,length(temp))];
    sample(j).freq_init(isnan(sample(j).freq_init))=0;
    select_freq=[select_freq,sample(j).freq_init(temp)];
    sample(j).rep_frac(isnan(sample(j).rep_frac))=0;
    select_rep=[select_rep,sample(j).rep_frac(temp)];
    end
    for q=1:length(select)
    [valuecorrmatrix(q,:),lags]=xcorr(sample(position(q)).fiber(select(q),:),sample(position(q)).fiber(select(q),:),'unbiased');
    end
    selectall{k}=select;
    positionall{k}=position;
    freqall{k}=select_freq;
    repall{k}=select_rep;
    valuecorrall{k}=valuecorrmatrix;
    k=k+1;
end

rf=cellfun(@mean,repall);
rferr=cellfun(@std,repall)./sqrt(cellfun(@length,repall));
inf=cellfun(@mean,freqall);
inferr=cellfun(@std,freqall)./sqrt(cellfun(@length,freqall));


for q=1:length(valuecorrall)
    if ~isempty(valuecorrall{q})
   funcorr(q,:)=mean(valuecorrall{q});
   funcorrerr(q,:)=std(valuecorrall{q})/sqrt(length(valuecorrall{q}(:,1)));
    else
        funcorr(q,:)=0;
   funcorrerr(q,:)=0;
    end
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

figure;
hold on
ylabel('I(t)');
xlabel('f(t)');
errorbar(rf,inf,inferr);
