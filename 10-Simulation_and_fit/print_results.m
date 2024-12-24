%I take the data from sample, I calculate the freq.of initiation and the
%correlation function in the indicated interval
%Then I fit the data with a genetic algorithm
%The speed is in kb/min, the frquency of initiation in 1/(kb*min)

%no automatic integer setting
clear all;
close all;

num_of_fits=100;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('results.mat');

for i=1:num_of_fits
  
vel_fit=results(i,2);
ceff=results(i,4);
a=results(i,3);

t(i,:)=(-(log(1-rf)*(ceff+1)*(ceff+2))/(2*a*vel_fit)).^(1/(ceff+2));
infth(i,:)=a.*t(i,:).^ceff;
tottime(i)=t(end);

Fiber_length=sample(1).Fiber_length;

for k=1:length(funcorr(:,1))
max_length1=2*vel_fit*t(i,k);
    if max_length1>Fiber_length-1
        max_length1=Fiber_length-1;
    end
x1=0:max_length1;
y1=funcorr(k,x1+Fiber_length);
phi1=exp(-2*vel_fit*a*t(i,k).^(ceff+2)/((ceff+1)*(ceff+2)));
comp1=t(i,k).*(1 - x1./max_length1);
xfit_unbiased1=1-2*phi1+(phi1^2).*exp(2*vel_fit*a.*comp1.^(ceff+2)/((ceff+1)*(ceff+2)));

allxfit_unbiased{k}{i}=xfit_unbiased1;
posi{k}{i}=x1;


end


end



scrsz = get(0,'ScreenSize');
figure('Position',[10 10 scrsz(3)*5/10 scrsz(4)*5/10])
subplot('Position',[0.15 0.15 0.35 0.8])
hold on
axis([0 ((rf(end)-rf(end-1))/2)+rf(end) 0 0.05])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
errorbar(rf,inf,inferr,'LineStyle','none','Marker','o','MarkerSize',2.5);
errorbar(rf,mean(infth),std(infth)/sqrt(num_of_fits),'r','linewidth',1);


subplot('Position',[0.6 0.15 0.35 0.8])
hold on
axis([0 40 0 1])
x=-Fiber_length+1:Fiber_length-1;
ylabel('C(r,f)','FontSize',12,'FontName','Arial');
xlabel('r (kb)','FontSize',12,'FontName','Arial');
for i=1:length(funcorr(:,1))
H=errorbar(x,funcorr(i,:),funcorrerr(i,:));
errorbar_tick(H,1,'units');
end

for k=1:length(funcorr(:,1))
interval=min(cellfun(@length,posi{k}));
a=cellfun(@(x) x(1:interval),allxfit_unbiased{k},'UniformOutput',false);
d=vertcat(a{:});
H1=errorbar(0:interval-1,mean(d),std(d)/sqrt(num_of_fits),'r','linewidth',1);
errorbar_tick(H,1,'units');
end

legend([H,H1],'Simulated data','Fit','FontSize',10);

im_paper1(['Result_fit_simulation'],6,3.2);

fprintf('Results fit simulation \n')
fprintf('\n')
fprintf('Avereage chi2 : %.2f \n',mean(results(:,1)))
fprintf('\n')
fprintf('Avereage parameters over the %i fits: \n',num_of_fits)
fprintf('v = %.4f; I0 = %.4f; alpha = %.4f \n',mean(results(:,2)),mean(results(:,3)),mean(results(:,4)))
fprintf('Average time to replicate a fiber from 0%% to %.2f%%: %.0f minutes \n',((rf(end)-rf(end-1))/2)+rf(end),mean(tottime))
fprintf('\n')


%close all

