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

load(['../5-Fit_two_processes/' sample_path '/results_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/rf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/minfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/maxfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorrerr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inferr_' sample_path '.mat']);


for i=1:num_of_fits
  
% The parameters in the function "fit_with_twoprocesses.m" are stored in a way 
% that the first three parameters results(i,2:4) correspond to the process
% with the higher speed; so in the code that follow the process labelled
% with 1 corresponds to the process with higher fork speed

vel=results(i,2);
ceff=results(i,4);
a=results(i,3);
vel2=results(i,5);
ceff2=results(i,7);
a2=results(i,6);
teta=results(i,8:end);


t(i,:)=(-(log(1-rf)*(ceff+1)*(ceff+2))/(2*a*vel)).^(1/(ceff+2));
t2(i,:)=(-(log(1-rf)*(ceff2+1)*(ceff2+2))/(2*a2*vel2)).^(1/(ceff2+2));
infth1(i,:)=a.*t(i,:).^ceff;
infth2(i,:)=a2.*t2(i,:).^ceff2;
infth(i,:)=teta.*infth1(i,:)+(1-teta).*infth2(i,:);
ttot1(i)=t(end); %number of steps
ttot2(i)=t2(end); %number of steps


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

max_length2=2*vel2*t2(i,k);
if max_length2>Fiber_length-1
    max_length2=Fiber_length-1;
end
x2=0:max_length2;
y2=funcorr(k,x2+Fiber_length);
phi2=exp(-2*vel2*a2*t2(i,k).^(ceff2+2)/((ceff2+1)*(ceff2+2)));
comp2=t2(i,k).*(1 - x2./max_length2);
xfit_unbiased2=1-2*phi2+(phi2^2).*exp(2*vel2*a2.*comp2.^(ceff2+2)/((ceff2+1)*(ceff2+2)));

if max_length1>max_length2
    if max_length1>Fiber_length-1
    max_length1=Fiber_length-1;
    end
    posi{k}{i}=x1;
    y=y1;
    xfit_unbiased2(end:floor(max_length1)+1)=xfit_unbiased2(end);
    xfit_unbiased{k}{i}=teta(k)*xfit_unbiased1+(1-teta(k))*xfit_unbiased2;
else
    if max_length2>Fiber_length-1
    max_length2=Fiber_length-1;
    end
    posi{k}{i}=x2;
    y=y2;
    xfit_unbiased1(end:floor(max_length2)+1)=xfit_unbiased1(end);
    xfit_unbiased{k}{i}=teta(k)*xfit_unbiased1+(1-teta(k))*xfit_unbiased2;
end


allxfit_unbiased1{k}{i}=xfit_unbiased1;
allxfit_unbiased2{k}{i}=xfit_unbiased2;
% xfit_biased=xfit_unbiased.*(maxfiberlength-x)./maxfiberlength;

end



end

scrsz = get(0,'ScreenSize');
figure('Position',[10 10 scrsz(3)*6/10 scrsz(4)*5/10])
subplot('Position',[0.10 0.15 0.23 0.8])
hold on
axis([0 ((rf(end)-rf(end-1))/2)+rf(end) 0 Inf])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
errorbar(rf,inf,inferr,'LineStyle','none','Marker','o','MarkerSize',2.5);
if num_of_fits==1
    plot(rf,infth,'r','linewidth',1);
    plot(rf,infth1,'g','linewidth',1);
    plot(rf,infth2,'k','linewidth',1);
else
    errorbar(rf,mean(infth),std(infth)/sqrt(num_of_fits),'r','linewidth',1);
    errorbar(rf,mean(infth1),std(infth2)/sqrt(num_of_fits),'g','linewidth',1);
    errorbar(rf,mean(infth2),std(infth2)/sqrt(num_of_fits),'k','linewidth',1);
end



subplot('Position',[0.40 0.15 0.28 0.8])
hold on
axis([0 minfiberlength 0 ((rf(end)-rf(end-1))/2)+rf(end)])
x=-maxfiberlength+1:maxfiberlength-1;
ylabel('C(r,f)','FontSize',12,'FontName','Arial');
xlabel('r (kb)','FontSize',12,'FontName','Arial');
for i=1:length(funcorr(:,1))
H=errorbar(x,funcorr(i,:),funcorrerr(i,:));
errorbar_tick(H,1,'units');
end

for k=1:length(funcorr(:,1))
interval=min(cellfun(@length,posi{k}));
a=cellfun(@(x) x(1:interval),xfit_unbiased{k},'UniformOutput',false);
b=cellfun(@(x) x(1:interval),allxfit_unbiased1{k},'UniformOutput',false);
c=cellfun(@(x) x(1:interval),allxfit_unbiased2{k},'UniformOutput',false);
d=vertcat(a{:});
e=vertcat(b{:});
f=vertcat(c{:});
if num_of_fits==1
    H1=plot(0:interval-1,d,'r','linewidth',1);
    H2=plot(0:interval-1,e,'g','linewidth',1);
    H3=plot(0:interval-1,f,'k','linewidth',1);
else
    H1=errorbar(0:interval-1,mean(d),std(d)/sqrt(num_of_fits),'r','linewidth',1);
    H2=errorbar(0:interval-1,mean(e),std(e)/sqrt(num_of_fits),'g','linewidth',1);
    H3=errorbar(0:interval-1,mean(f),std(f)/sqrt(num_of_fits),'k','linewidth',1);
end
errorbar_tick(H,1,'units');
end

legend([H,H1,H2,H3],'Data','Fit','Process 1','Process 2');

subplot('Position',[0.76 0.15 0.23 0.8])
hold on
axis([0 ((rf(end)-rf(end-1))/2)+rf(end) 0 Inf])
ylabel('\theta','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
if num_of_fits==1
    plot(rf,results(:,8:end),'r','linewidth',1);
else
    errorbar(rf,mean(results(:,8:end)),std(results(:,8:end))/sqrt(num_of_fits),'r','linewidth',1);
end


im_paper1([sample_path '/Result_fit_' sample_path ],8,3.2);

fprintf('Results %s \n',sample_path)
fprintf('\n')
fprintf('Avereage chi2 : %.2f \n',mean(results(:,1)))
fprintf('\n')
fprintf('Avereage parameters over the %i fits: \n',num_of_fits)
fprintf('Process 1: v = %.4f; I0 = %.4f; alpha = %.4f \n',mean(results(:,2)),mean(results(:,3)),mean(results(:,4)))
fprintf('Process 2: v = %.4f; I0 = %.4f; alpha = %.4f \n',mean(results(:,5)),mean(results(:,6)),mean(results(:,7)))
fprintf('Average time to replicate a fiber from 0%% to %.2f%% for process1: %.0f minutes \n',((rf(end)-rf(end-1))/2)+rf(end),mean(ttot1))
fprintf('Average time to replicate a fiber from 0%% to %.2f%% for process2: %.0f minutes \n',((rf(end)-rf(end-1))/2)+rf(end),mean(ttot2))
fprintf('\n')
fprintf('Avereage theta parameters: \n')
fprintf('%.4f ',mean(results(:,8:end)))
fprintf('\n')

close all