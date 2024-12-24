%_______________________ANALYSIS OF THE RESULTS___________________________
%
% With this program I represent the results of the fits with the genetic
% algorithm
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%

sample_path='output_demo';
num_of_fits=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['../4-Fit_single_process/' sample_path '/results_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/rf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/minfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/maxfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorrerr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inferr_' sample_path '.mat']);


for i=1:num_of_fits
  
vel=results(i,2);
ceff=results(i,4);
a=results(i,3);

t(i,:)=(-(log(1-rf)*(ceff+1)*(ceff+2))/(2*a*vel)).^(1/(ceff+2));
infth(i,:)=a.*t(i,:).^ceff;
tottime(i)=t(end);

Fiber_length=maxfiberlength;

for k=1:length(funcorr(:,1))
max_length1=2*vel*t(i,k);
    if max_length1>Fiber_length-1
        max_length1=Fiber_length-1;
    end
x1=0:max_length1;
y1=funcorr(k,x1+Fiber_length);
phi1=exp(-2*vel*a*t(i,k).^(ceff+2)/((ceff+1)*(ceff+2)));
comp1=t(i,k).*(1 - x1./max_length1);
xfit_unbiased1=1-2*phi1+(phi1^2).*exp(2*vel*a.*comp1.^(ceff+2)/((ceff+1)*(ceff+2)));

allxfit_unbiased{k}{i}=xfit_unbiased1;
posi{k}{i}=x1;


end


end



scrsz = get(0,'ScreenSize');
figure('Position',[10 10 scrsz(3)*5/10 scrsz(4)*5/10])
subplot('Position',[0.15 0.15 0.35 0.8])
hold on
axis([0 ((rf(end)-rf(end-1))/2)+rf(end) 0 Inf])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
errorbar(rf,inf,inferr,'LineStyle','none','Marker','o','MarkerSize',2.5);
if num_of_fits==1
    plot(rf,infth,'r','linewidth',1);
else
    errorbar(rf,mean(infth),std(infth)/sqrt(num_of_fits),'r','linewidth',1);
end

subplot('Position',[0.6 0.15 0.35 0.8])
hold on
axis([0 40 0 1])
x=-maxfiberlength+1:maxfiberlength-1;
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
if num_of_fits==1
    H1=plot(0:interval-1,d,'r','linewidth',1);
else
    H1=errorbar(0:interval-1,mean(d),std(d)/sqrt(num_of_fits),'r','linewidth',1);
end
errorbar_tick(H,1,'units');
end

legend([H,H1],'Experiment','Fit');

im_paper1([sample_path '/Result_fit_oneprocess_' sample_path ],6,3.2);

fprintf('Results %s \n',sample_path)
fprintf('\n')
fprintf('Avereage chi2 : %.2f \n',mean(results(:,1)))
fprintf('\n')
fprintf('Avereage parameters over the %i fits: \n',num_of_fits)
fprintf('v = %.4f; I0 = %.10f; alpha = %.4f \n',mean(results(:,2)),mean(results(:,3)),mean(results(:,4)))
fprintf('Average time to replicate a fiber from 0%% to %.2f%%: %.0f minutes \n',((rf(end)-rf(end-1))/2)+rf(end),mean(tottime))
fprintf('\n')


close all