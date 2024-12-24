%_______________________ANALYSIS OF THE RESULTS___________________________
%
% With this program I compare the results between control and treatment
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%

sample_path_tosave='output_demo';
sample_path1='output_demo';
sample_path2='output_demo_treatment';
num_of_fits=100;

%%%%%%%%%%%%%%%%%%%%%%%Comparison%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['../2-Plot_autocorrelation/' sample_path1 '/rf_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/minfiberlength_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/maxfiberlength_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/funcorr_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/funcorrerr_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/inf_' sample_path1 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path1 '/inferr_' sample_path1 '.mat']);
load(['../5-Fit_two_processes/' sample_path1 '/resultsFINAL_' sample_path1 '.mat']);

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

%average1=mean(results)
%error1=std(results)./10
%time1=[mean(ttot1) std(ttot1)/10]
%time2=[mean(ttot2) std(ttot2)/10]


scrsz = get(0,'ScreenSize');
s=figure('Position',[10 10 scrsz(3)*6/10 scrsz(4)*5/10]);
s1=subplot('Position',[0.10 0.15 0.23 0.8]);
hold on
axis([0 0.75 0 Inf])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
plot(rf,mean(infth1),'g','linewidth',1);
plot(rf,mean(infth2),'k','linewidth',1);
infth2_forratio(1,:)=mean(infth2);
infth1_forratio(1,:)=mean(infth1);


s2=subplot('Position',[0.40 0.15 0.28 0.8]);
hold on
axis([0 minfiberlength 0 0.75])
x=-maxfiberlength+1:maxfiberlength-1;
ylabel('C(r,f)','FontSize',12,'FontName','Arial');
xlabel('r (kb)','FontSize',12,'FontName','Arial');

for k=1:length(funcorr(:,1))
% H1=plot(mean(x{k}),xfit_unbiased,'r','linewidth',1);
% H2=plot(x,xfit_unbiased1,'g','linewidth',1);
% H3=plot(x,xfit_unbiased2,'k','linewidth',1);
interval=min(cellfun(@length,posi{k}));
a=cellfun(@(x) x(1:interval),xfit_unbiased{k},'UniformOutput',false);
b=cellfun(@(x) x(1:interval),allxfit_unbiased1{k},'UniformOutput',false);
c=cellfun(@(x) x(1:interval),allxfit_unbiased2{k},'UniformOutput',false);
d=vertcat(a{:});
e=vertcat(b{:});
f=vertcat(c{:});
H1=plot(0:interval-1,mean(e),'g','linewidth',1);
H2=plot(0:interval-1,mean(f),'k','linewidth',1);
% errorbar_tick(H1,1,'units');
% errorbar_tick(H2,1,'units');
end


s3=subplot('Position',[0.76 0.15 0.23 0.8]);
hold on
axis([0 0.75 0 Inf])
ylabel('\theta','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
plot(rf,mean(results(:,8:end)),'r','linewidth',1);

clear t t2 ttot1 ttot2 infth1 infth2 infth teta posi xfit_unbiased allxfit_unbiased1 allxfit_unbiased2
load(['../2-Plot_autocorrelation/' sample_path2 '/rf_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/minfiberlength_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/maxfiberlength_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/funcorr_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/funcorrerr_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/inf_' sample_path2 '.mat']);
load(['../2-Plot_autocorrelation/' sample_path2 '/inferr_' sample_path2 '.mat']);
load(['../5-Fit_two_processes/' sample_path2 '/resultsFINAL_' sample_path2 '.mat']);

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
figure(s);
subplot(s1);
hold on
axis([0 0.75 0 Inf])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
plot(rf,mean(infth1),'-og','linewidth',1,'markersize',2);
plot(rf,mean(infth2),'-ok','linewidth',1,'markersize',2);
infth2_forratio(2,:)=mean(infth2);
infth1_forratio(2,:)=mean(infth1);
%Ratio I(f)_treatment/I(f)_control
fprintf('Ratio I(f)_treat/I(f)_control process fast \n')
fprintf('%.4f ',infth1_forratio(2,:)./infth1_forratio(1,:))
fprintf('\n')
fprintf('Ratio I(f)_treat/I(f)_control process slow \n')
fprintf('%.4f ',infth2_forratio(2,:)./infth2_forratio(1,:))
fprintf('\n')

subplot(s2);
hold on
axis([0 50 0 0.75])
x=-maxfiberlength+1:maxfiberlength-1;
ylabel('C(r,f)','FontSize',12,'FontName','Arial');
xlabel('r (kb)','FontSize',12,'FontName','Arial');


for k=1:length(funcorr(:,1))
interval=min(cellfun(@length,posi{k}));
a=cellfun(@(x) x(1:interval),xfit_unbiased{k},'UniformOutput',false);
b=cellfun(@(x) x(1:interval),allxfit_unbiased1{k},'UniformOutput',false);
c=cellfun(@(x) x(1:interval),allxfit_unbiased2{k},'UniformOutput',false);
d=vertcat(a{:});
e=vertcat(b{:});
f=vertcat(c{:});
H3=plot(0:interval-1,mean(e),'-og','linewidth',1,'markersize',2);
H4=plot(0:interval-1,mean(f),'-ok','linewidth',1,'markersize',2);
% errorbar_tick(H3,1,'units');
% errorbar_tick(H4,1,'units');
end

legend([H1,H2,H3,H4],'Process 1 Ctrl','Process 2 Ctrl','Process 1 Depl','Process 2 Depl','FontSize',9);

subplot(s3);
hold on
axis([0 0.75 0 Inf])
ylabel('\theta','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
plot(rf,mean(results(:,8:end)),'-or','linewidth',1,'markersize',2);
legend({'Ctrl' 'Depl'},'FontSize',10);

im_paper1([sample_path_tosave '/Result_fit_' sample_path_tosave ],8,3.2);
