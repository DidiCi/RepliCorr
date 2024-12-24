
%___________________________CLUSTERING AND COMPARISON____________________
%
% With this program I analyse more in detail the results of the fit and I
% classify the fibers depending on the result of the fit
%________________________________________________________________________

clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%

sample_path='output_demo';
num_of_fits=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['../5-Fit_two_processes/' sample_path '/results_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/sample_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/rf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/minfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/maxfiberlength_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/funcorrerr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inf_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/inferr_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/positionall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/freqall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/repall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/valuecorrall_' sample_path '.mat'])
load(['../2-Plot_autocorrelation/' sample_path '/eyesall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/fibersall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/etedall_' sample_path '.mat']);
load(['../3-PCA_and_clustering/' sample_path '/T_' sample_path '.mat']);


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

% The correlation profiles are calculated by using the parameters 
% of the fit

Fiber_length=maxfiberlength;
for k=1:length(funcorr(:,1))
    max_length1=2*vel*t(i,k);
    if max_length1>Fiber_length-1
        max_length1=Fiber_length-1;
    end
    x1=0:max_length1;
    phi1=exp(-2*vel*a*t(i,k).^(ceff+2)/((ceff+1)*(ceff+2)));
    comp1=t(i,k).*(1 - x1./max_length1);
    xfit_unbiased1=1-2*phi1+(phi1^2).*exp(2*vel*a.*comp1.^(ceff+2)/((ceff+1)*(ceff+2)));

    max_length2=2*vel2*t2(i,k);
    if max_length2>Fiber_length-1
        max_length2=Fiber_length-1;
    end
    x2=0:max_length2;
    phi2=exp(-2*vel2*a2*t2(i,k).^(ceff2+2)/((ceff2+1)*(ceff2+2)));
    comp2=t2(i,k).*(1 - x2./max_length2);
    xfit_unbiased2=1-2*phi2+(phi2^2).*exp(2*vel2*a2.*comp2.^(ceff2+2)/((ceff2+1)*(ceff2+2)));

    % "max_length1" and "max_length2" are compared and the largest of the
    % two is used as limit for the x of the theoretical correlation
    % function. The other function is completed by repeating the last value
    % of the function for the missing length.
   
    if max_length1>max_length2 
        posi{k}{i}=x1;
        xfit_unbiased2(end:floor(max_length1)+1)=xfit_unbiased2(end);
        xfit_unbiased{k}{i}=teta(k)*xfit_unbiased1+(1-teta(k))*xfit_unbiased2;
    else
        posi{k}{i}=x2;
        xfit_unbiased1(end:floor(max_length2)+1)=xfit_unbiased1(end);
        xfit_unbiased{k}{i}=teta(k)*xfit_unbiased1+(1-teta(k))*xfit_unbiased2;
    end

    % k is used to label the different bins of replicated fraction.
    % i is used to label the different fits.
    % For each bin of replicated fraction, the limit of the x axis for 
    % the theoretical correlation function is chosen as the minimum
    % of the different "max_length" of the bin
    
    allxfit_unbiased1{k}{i}=xfit_unbiased1;
    allxfit_unbiased2{k}{i}=xfit_unbiased2;
    x1all{k}=0:min(cellfun(@length,posi{k}))-1;
    
end

end


% Loop on the bins of replicated fraction
for j=1:length(funcorr(:,1))
    
    % The parameters in the function "fit_with_twoprocesses.m" are stored in a way 
    % that the first three parameters results(i,2:4) correspond to the process
    % with the higher speed; so in the code that follow the process labelled
    % with 1 corresponds to the process with higher fork speed

    x1=x1all{j};
    aa=cellfun(@(x) x(x1+1),xfit_unbiased{j},'UniformOutput',false);
    bb=cellfun(@(x) x(x1+1),allxfit_unbiased1{j},'UniformOutput',false);
    cc=cellfun(@(x) x(x1+1),allxfit_unbiased2{j},'UniformOutput',false);
    if num_of_fits==1
        xfit_unb=vertcat(aa{:});
        xfit_unb1=vertcat(bb{:}); 
        xfit_unb2=vertcat(cc{:}); 
    else
        xfit_unb=mean(vertcat(aa{:})); % Average total correlation function
        xfit_unb1=mean(vertcat(bb{:})); % Average total correlation function for process 1
        xfit_unb2=mean(vertcat(cc{:})); % Average total correlation function for process 2
    end


    % Correlation between model and experiments; X, Y and Y2 contain
    % the correlation profiles
    
    for q=1:length(valuecorrall{j}(:,1))
        X{j}(q,:)=valuecorrall{j}(q,x1+Fiber_length)./valuecorrall{j}(q,Fiber_length);
    end

    Y{j}=xfit_unb1./xfit_unb1(1);
    Y2{j}=xfit_unb2./xfit_unb2(1);
    
    % Calculate pairwise distance between two sets of observations
    
    D{j} = pdist2(X{j},Y{j},'correlation'); %small value means correlated
    D2{j} = pdist2(X{j},Y2{j},'correlation');

    figure(9); subplot(1,7,j); 
    hold on; 
    axis([0 Inf 0 Inf])
    title(['Bin' num2str(j) ' (f=' num2str(round(rf(j)*100)/100) ')'],'FontSize',14,'FontName','Arial');
    scatter(D{j}(T{j}==1),D2{j}(T{j}==1),'r','filled'); 
    scatter(D{j}(T{j}==2),D2{j}(T{j}==2),'k','filled');
    ylabel('Distance slow process','FontSize',14,'FontName','Arial');
    xlabel('Distance fast process','FontSize',14,'FontName','Arial');
    if j==1
        legend('Cluster 1','Cluster 2');
    end
    plot([0 1],[0 1]);

    figure(10); subplot(1,7,j); 
    hold on; 
    color_points=[0.0 0.4 1.0];
    color_line=[1.0 0.2 0.2];
    axis([0 Inf 0 Inf])
    title(['Bin' num2str(j) ' (f=' num2str(round(rf(j)*100)/100) ')'],'FontSize',14,'FontName','Arial');
    scatter(D{j}(D{j}<D2{j}),D2{j}(D{j}<D2{j}),15,color_points,'filled'); 
    scatter(D{j}(D{j}>D2{j}),D2{j}(D{j}>D2{j}),15,color_points,'filled');
    ylabel('Distance slow process','FontSize',14,'FontName','Arial');
    xlabel('Distance fast process','FontSize',14,'FontName','Arial');
    plot([0 1],[0 1],'color',color_line,'LineWidth',1.5);
    points_transition_repfiber{j}.x=D{j};
    points_transition_repfiber{j}.y=D2{j};
    save([sample_path '/points_transition_repfiber_' sample_path '.mat'], 'points_transition_repfiber');
    
    figure(11);
    hold on; 
    scatter(rf(j),sum(D{j}<D2{j})/length(D{j}),15,color_points); 
    axis([0 Inf 0 1])
    ylabel('Fraction fast process','FontSize',14,'FontName','Arial');
    xlabel('Replicated fraction of fiber','FontSize',14,'FontName','Arial');
    title('Fraction of fibers following fast process by replicated fraction of fiber')
    
      
    %I use the distance in D2 because in the other category same fibers below
    %and abore the curve are considered far
    %I set for each fiber:
    % - category 1 for fibers following fast process
    % - category 2 for fibers following slow process
    cat1{j}=find(D{j}<D2{j});
    cat2{j}=find(D{j}>=D2{j});
    num1(j)=length(cat1{j});
    num2(j)=length(cat2{j});
    %I plot the fiber in the two categories with two different colors
    figure(5); subplot(1,7,j); hold on; 
    ylabel('C(r,f)/C(0,f)','FontSize',12,'FontName','Arial');
    xlabel('r (kb)','FontSize',12,'FontName','Arial');
    title(['Bin' num2str(j) ' (f=' num2str(round(rf(j)*100)/100) ')'],'FontSize',12,'FontName','Arial');
    axis([0 length(Y{j})+1 -Inf Inf])
    % plot(Y{j},'k'); plot(Y2{j},'k'); 
    for q=1:length(cat1{j})
        plot(X{j}(cat1{j}(q),:),'g');
    end
    for q=1:length(cat2{j})
        plot(X{j}(cat2{j}(q),:),'r')
    end

    % I plot the averages for the two categories
    if isempty(cat1{j})
        firstcat=zeros(1,length(Y{j}));
        firstcaterr=zeros(1,length(Y{j}));
        secondcat=mean(valuecorrall{j}(cat2{j},x1+Fiber_length));
        secondcaterr=std(valuecorrall{j}(cat2{j},x1+Fiber_length))/sqrt(length(cat2{j}));
    elseif length(cat1{j})==length(D{j})
        firstcat=mean(valuecorrall{j}(cat1{j},x1+Fiber_length));
        firstcaterr=std(valuecorrall{j}(cat1{j},x1+Fiber_length))/sqrt(length(cat1{j}));
        secondcat=zeros(1,length(Y{j}));
        secondcaterr=zeros(1,length(Y{j}));
    elseif length(cat1{j})==1
        firstcat=valuecorrall{j}(cat1{j},x1+Fiber_length);
        firstcaterr=zeros(1,length(Y{j}));
        secondcat=mean(valuecorrall{j}(cat2{j},x1+Fiber_length));
        secondcaterr=std(valuecorrall{j}(cat2{j},x1+Fiber_length))/sqrt(length(cat2{j}));
    elseif length(cat2{j})==1
        firstcat=mean(valuecorrall{j}(cat1{j},x1+Fiber_length));
        firstcaterr=std(valuecorrall{j}(cat1{j},x1+Fiber_length))/sqrt(length(cat1{j}));
        secondcat=valuecorrall{j}(cat2{j},x1+Fiber_length);   
        secondcaterr=zeros(1,length(Y{j}));
    else
        firstcat=mean(valuecorrall{j}(cat1{j},x1+Fiber_length));
        firstcaterr=std(valuecorrall{j}(cat1{j},x1+Fiber_length))/sqrt(length(cat1{j}));
        secondcat=mean(valuecorrall{j}(cat2{j},x1+Fiber_length));
        secondcaterr=std(valuecorrall{j}(cat2{j},x1+Fiber_length))/sqrt(length(cat2{j}));
    end

    figure(6);
    subplot(1,8,j); 
    hold on; 
    title(['Bin' num2str(j) ' (f=' num2str(round(rf(j)*100)/100) ')'],'FontSize',12,'FontName','Arial');
    axis([0 length(firstcat)+1 -Inf Inf])
    x=0:length(firstcat)-1;
    ylabel('C(r,f)','FontSize',12,'FontName','Arial');
    xlabel('r (kb)','FontSize',12,'FontName','Arial');
    H=errorbar(x,firstcat,firstcaterr,'r');
    errorbar_tick(H,1,'units');
    H1=errorbar(x,secondcat,secondcaterr,'b');
    errorbar_tick(H,1,'units');
    plot(x,xfit_unb1,'g','linewidth',1); 
    plot(x,xfit_unb2,'k','linewidth',1); 
    hold off
    if j==1
        legend('Category 1','Category 2','Model 1','Model 2');
    end


    %I calculate averaged frquency of initiation and error
    freq1(1,j)=mean(freqall{j}(cat1{j}));
    freq1(2,j)=std(freqall{j}(cat1{j}))/sqrt(length(cat1{j}));
    %freq1(2,j)=max_error(freqall{j}(cat1{j}));
    freq2(1,j)=mean(freqall{j}(cat2{j}));
    freq2(2,j)=std(freqall{j}(cat2{j}))/sqrt(length(cat2{j}));
    %freq2(2,j)=max_error(freqall{j}(cat2{j}));
    
    %I calculate averaged replicated fraction and error
    rep1(1,j)=mean(repall{j}(cat1{j}));
    rep1(2,j)=std(repall{j}(cat1{j}))/sqrt(length(cat1{j}));
    rep2(1,j)=mean(repall{j}(cat2{j}));
    rep2(2,j)=std(repall{j}(cat2{j}))/sqrt(length(cat2{j}));
    
    %I calculate averaged eyes length and error
    eyes1(1,j)=mean([eyesall{j}{cat1{j}}]);
    eyes1(2,j)=std([eyesall{j}{cat1{j}}])/sqrt(length([eyesall{j}{cat1{j}}]));
    %eyes1(2,j)=max_error([eyesall{j}{cat1{j}}]);
    eyes2(1,j)=mean([eyesall{j}{cat2{j}}]);
    eyes2(2,j)=std([eyesall{j}{cat2{j}}])/sqrt(length([eyesall{j}{cat2{j}}]));
    %eyes2(2,j)=max_error([eyesall{j}{cat2{j}}]);

end

%I plot the fraction of fibers
figure(6);
subplot(1,8,8);
hold on
ylabel('Fraction of fibers');
xlabel('Bin');
plot(1-teta,'k','linewidth',1);
plot(num1./(num1+num2),'r','linewidth',1);
plot(teta,'g','linewidth',1);
plot(num2./(num1+num2),'b','linewidth',1);


figure(7);
subplot(1,2,1)
hold on;  
axis([0 0.75 -Inf Inf])
ylabel('I(f) (1/(kb*min))','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
errorbar(rep1(1,:),freq1(1,:),freq1(2,:),'r','LineStyle','none','Marker','o','MarkerSize',2.5);
errorbar(rep2(1,:),freq2(1,:),freq2(2,:),'b','LineStyle','none','Marker','o','MarkerSize',2.5);
plot(rf,infth1,'g','linewidth',1);
plot(rf,infth2,'k','linewidth',1);
subplot(1,2,2)
hold on;  
axis([0 0.75 -Inf Inf])
ylabel('Mean eye lengths (kb)','FontSize',12,'FontName','Arial');
xlabel('f','FontSize',12,'FontName','Arial');
errorbar(rep1(1,:),eyes1(1,:),eyes1(2,:),'r','LineStyle','none','Marker','o','MarkerSize',2.5);
errorbar(rep2(1,:),eyes2(1,:),eyes2(2,:),'b','LineStyle','none','Marker','o','MarkerSize',2.5);
plot(rf,2*vel*t,'g','linewidth',1);
plot(rf,2*vel2*t2,'k','linewidth',1);

for num_bin=1:length(cat1)
    
    figure;
    subplot(2,1,1)
    hold on; 
    nice_green=[0.2,0.8,0.2];
    nice_red=[1,0.0,0.2];
    num_fibers_toprint=10;
    Bin=num_bin;
    axis([0 Inf 0 num_fibers_toprint*2+3])
    title(['Fast process - Replicated fraction ' num2str(rep1(1,Bin))]);
    %for i=1:length(cat1{Bin1})
    for i=1:num_fibers_toprint
        try
            red=find(fibersall{Bin}{cat1{Bin}(i)}==1)';
            green=find(fibersall{Bin}{cat1{Bin}(i)}==0)';
            scatter(red,2*i.*ones(1,length(red)),10,nice_red,'s','filled');
            scatter(green,2*i.*ones(1,length(green)),10,nice_green,'s','filled');
        catch
            fprintf('Not enough fibers in Bin: %i Process: Fast\n', Bin);
            break;  
        end
    end
    xlabel('Fiber length (kb)','FontSize',12,'FontName','Arial');
    set(gca,'YTickLabel',[]);
    set(gca,'YColor','w');
    
    subplot(2,1,2)
    hold on;
    axis([0 Inf 0 num_fibers_toprint*2+3])
    title(['Slow process - Replicated fraction ' num2str(rep1(1,Bin))]);
    for i=1:num_fibers_toprint
        try
            red=find(fibersall{Bin}{cat2{Bin}(i)}==1)';
            green=find(fibersall{Bin}{cat2{Bin}(i)}==0)';
            scatter(red,2*i.*ones(1,length(red)),10,nice_red,'s','filled');
            scatter(green,2*i.*ones(1,length(green)),10,nice_green,'s','filled');
        catch
            fprintf('Not enough fibers in Bin: %i Process: Slow\n', Bin);
            break;  
        end
    end
    xlabel('Fiber length (kb)','FontSize',12,'FontName','Arial');
    set(gca,'YTickLabel',[]);
    set(gca,'YColor','w');
    im_paper1([sample_path '/Example_fibers_' sample_path '_Bin' num2str(rep1(1,Bin))],5,6);
    
    xlswrite([sample_path '/Fibers_rawimages_fastprocess.xls'],{'Raw fibers by bin of replicated fraction and process',['Bin ' num2str(rep1(1,Bin))]},num_bin,'A1:B1');
    xlswrite([sample_path '/Fibers_rawimages_fastprocess.xls'],{'Experiment','File number in log.txt','Fiber id'},num_bin,'A2:B2');
    xlswrite([sample_path '/Fibers_rawimages_fastprocess.xls'],{sample(positionall{Bin}(cat1{Bin})).file}',num_bin,['A3:A' num2str(2+length({sample(positionall{Bin}(cat1{Bin})).file}))]);
    xlswrite([sample_path '/Fibers_rawimages_fastprocess.xls'],[ sample(positionall{Bin}(cat1{Bin})).position]',num_bin,['B3:B' num2str(2+length({sample(positionall{Bin}(cat1{Bin})).file}))]);
    xlswrite([sample_path '/Fibers_rawimages_fastprocess.xls'],{ sample(positionall{Bin}(cat1{Bin})).fiber_id}',num_bin,['C3:C' num2str(2+length({sample(positionall{Bin}(cat1{Bin})).file}))]);

    xlswrite([sample_path '/Fibers_rawimages_slowprocess.xls'],{'Raw fibers by bin of replicated fraction and process',['Bin ' num2str(rep1(1,Bin))]},num_bin,'A1:B1');
    xlswrite([sample_path '/Fibers_rawimages_slowprocess.xls'],{'Experiment','File number in log.txt','Fiber id'},num_bin,'A2:B2');
    xlswrite([sample_path '/Fibers_rawimages_slowprocess.xls'],{sample(positionall{Bin}(cat2{Bin})).file}',num_bin,['A3:A' num2str(2+length({sample(positionall{Bin}(cat2{Bin})).file}))]);
    xlswrite([sample_path '/Fibers_rawimages_slowprocess.xls'],[ sample(positionall{Bin}(cat2{Bin})).position]',num_bin,['B3:B' num2str(2+length({sample(positionall{Bin}(cat2{Bin})).file}))]);
    xlswrite([sample_path '/Fibers_rawimages_slowprocess.xls'],{sample(positionall{Bin}(cat2{Bin})).fiber_id}',num_bin,['C3:C' num2str(2+length({sample(positionall{Bin}(cat2{Bin})).file}))]);

end

figure;
%Number of bins by which the progam divide the fibers according to
%the replicated fraction (I didn't use the fuction 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_binrep=8;
%Number of bins by which the progam divide the ETED, eye length and gap length
%distributions(I didn't use the fuction 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_bineyes=20;
%Max value by which the progam divide the ETED, eye length and gap length
%distributions
maxlength_bineyes=80;
centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);

histeted=zeros(1,length(centerseyes));
histeyes=zeros(1,length(centerseyes));
histgaps=zeros(1,length(centerseyes));

k=1;
%I extract the ETED, eye lengths and gap lengths from all the fibers of the
%given time point.
%There will be no 0 values (because are all lengths) so I don't need to add 
%the 0 values in the first bin like for the replicated fraction
eted1=[];
eted2=[];
for i=1:length(cat1)
eted1=[eted1,[etedall{i}{cat1{i}}]];
eted2=[eted2,[etedall{i}{cat2{i}}]];
end
for i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes
    %I calculate the number of ETED in the bin
    histeted1(k)=sum(eted1>i & eted1<=i+(maxlength_bineyes/(num_bineyes)));
    histeted2(k)=sum(eted2>i & eted2<=i+(maxlength_bineyes/(num_bineyes)));
    k=k+1;
end
%I normalize the distributions and I put to 0 the NaN
%values; a point is NaN if there are not origins or forks for that bin.
histeted1=histeted1./sum(histeted1);
histeted2=histeted2./sum(histeted2);
histeted1(isnan(histeted1))=0;
histeted2(isnan(histeted2))=0;

hold on
plot(centerseyes,histeted1,'b');
plot(centerseyes,histeted2,'r');
legend('Model1','Model2');



%theta faction by tot replicated fraction
%I set for each fiber:
% - category 1 for fibers following fast process
% - category 2 for fibers following slow process
% - category 0 for not replicated and completely replicated fibers

for j=1:length(funcorr(:,1))
    
    [sample(positionall{j}(cat1{j})).category]=deal(1);
    [sample(positionall{j}(cat2{j})).category]=deal(2);
    
    [sample(cellfun(@isempty, {sample.category})).category]=deal(0);
    Dassign = num2cell(D{j});
    D2assign = num2cell(D2{j});
    [sample(positionall{j}).D]=Dassign{:};
    [sample(positionall{j}).D2]=D2assign{:};
end

%One way to select fast or low fibers by total replicated fraction...
expfrac_tot=unique([sample.totrep_frac]);
for i=1:length(expfrac_tot)
    fibercat1_totrep(i)=sum((([sample.totrep_frac]==expfrac_tot(i)) & ([sample.category]==1)));
    fibercat2_totrep(i)=sum((([sample.totrep_frac]==expfrac_tot(i)) & ([sample.category]==2)));
end


figure;
hold on
axis([0 Inf 0 1])
scatter(expfrac_tot,fibercat1_totrep./(fibercat2_totrep+fibercat1_totrep),15,color_points)
ylabel('Fraction fast process','FontSize',14,'FontName','Arial');
xlabel('Total replicated fraction','FontSize',14,'FontName','Arial');
title('Fraction of fibers following fast process by total rep.frac.')
im_paper1([sample_path '/FractFiber_FastProcess_TotRep_' sample_path],5,3);


figure;
hold on
scatter(expfrac_tot,fibercat1_totrep,'r')
scatter(expfrac_tot,fibercat2_totrep,'b')
legend('num fast','num slow')
title('Number of fibers following processes by total rep.frac.')

figure;

for i=1:length(expfrac_tot)
hold on
axis([0 Inf 0 1])

fraction_process1=sum([sample([sample.totrep_frac]==expfrac_tot(i)).D] < [sample([sample.totrep_frac]==expfrac_tot(i)).D2])/length([sample([sample.totrep_frac]==expfrac_tot(i)).D]);
scatter(expfrac_tot(i),fraction_process1,40,color_points,'filled'); 

%scatter(D{j}(D{j}<D2{j}),D2{j}(D{j}<D2{j}),'g','filled'); 
%scatter(D{j}(D{j}>D2{j}),D2{j}(D{j}>D2{j}),'r','filled');
% scatter(D{j}(T{j}==3),D2{j}(T{j}==3),'g','filled');
% scatter(D{j}(T{j}==4),D2{j}(T{j}==4),'k','filled');
% scatter(D{j}(T{j}==5),D2{j}(T{j}==5),'y','filled');
ylabel('Fraction fast process','FontSize',14,'FontName','Arial');
xlabel('Total replicated fraction','FontSize',14,'FontName','Arial');
title('Fraction of fibers following fast process by total rep.frac.')
% if j==1
% legend('Cluster 1','Cluster 2');
% end
hold off
end

%Each plot correspond to a time point;
%colors by replicated fraction of the fiber
figure;

for i=1:length(expfrac_tot)
    subplot(1,length(expfrac_tot),i)
 hold on
axis([0 Inf 0 Inf])
title(['Bin' num2str(i) ' (tot_f=' num2str(round(expfrac_tot(i)*100)/100) ')'],'FontSize',14,'FontName','Arial');

points_early_D=[sample(([sample.fracfiber]>0) & ([sample.fracfiber]<=0.25) & ([sample.totrep_frac]==expfrac_tot(i))).D];
points_early_D2=[sample(([sample.fracfiber]>0) & ([sample.fracfiber]<=0.25) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
H1=scatter(points_early_D,points_early_D2,15,[0 0.4470 0.7410],'filled'); 
points_middle_D=[sample(([sample.fracfiber]>0.25) & ([sample.fracfiber]<=0.5) & ([sample.totrep_frac]==expfrac_tot(i))).D];
points_middle_D2=[sample(([sample.fracfiber]>0.25) & ([sample.fracfiber]<=0.5) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
H2=scatter(points_middle_D,points_middle_D2,15,'k','filled'); 
points_late_D=[sample(([sample.fracfiber]>0.5) & ([sample.fracfiber]<=0.75) & ([sample.totrep_frac]==expfrac_tot(i))).D];
points_late_D2=[sample(([sample.fracfiber]>0.5) & ([sample.fracfiber]<=0.75) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
H3=scatter(points_late_D,points_late_D2,15,[0.8500 0.3250 0.0980],'filled'); 

%scatter(D{j}(D{j}<D2{j}),D2{j}(D{j}<D2{j}),'g','filled'); 
%scatter(D{j}(D{j}>D2{j}),D2{j}(D{j}>D2{j}),'r','filled');
% scatter(D{j}(T{j}==3),D2{j}(T{j}==3),'g','filled');
% scatter(D{j}(T{j}==4),D2{j}(T{j}==4),'k','filled');
% scatter(D{j}(T{j}==5),D2{j}(T{j}==5),'y','filled');
ylabel('Distance slow process','FontSize',14,'FontName','Arial');
xlabel('Distance fast process','FontSize',14,'FontName','Arial');
if i==1
    legend('Low rep. (0-25%)','Middle rep. (25-50%)','High rep. (50-75%)');
end
plot([0 1],[0 1],'color',color_line,'LineWidth',1.5);
hold off
end
im_paper1([sample_path '/Transition_TotRep_byRepFiber' sample_path],20,3);

%Proportion of points in the previous figure
figure;
for i=1:length(expfrac_tot)
    hold on; 
    points_early_D=[sample(([sample.fracfiber]>0) & ([sample.fracfiber]<=0.25) & ([sample.totrep_frac]==expfrac_tot(i))).D];
    points_early_D2=[sample(([sample.fracfiber]>0) & ([sample.fracfiber]<=0.25) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
    H1=scatter(expfrac_tot(i),sum(points_early_D<points_early_D2)/length(points_early_D),'g','filled'); 
    points_middle_D=[sample(([sample.fracfiber]>0.25) & ([sample.fracfiber]<=0.5) & ([sample.totrep_frac]==expfrac_tot(i))).D];
    points_middle_D2=[sample(([sample.fracfiber]>0.25) & ([sample.fracfiber]<=0.5) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
    H2=scatter(expfrac_tot(i),sum(points_middle_D<points_middle_D2)/length(points_middle_D),'y','filled'); 
    points_late_D=[sample(([sample.fracfiber]>0.5) & ([sample.fracfiber]<=0.75) & ([sample.totrep_frac]==expfrac_tot(i))).D];
    points_late_D2=[sample(([sample.fracfiber]>0.5) & ([sample.fracfiber]<=0.75) & ([sample.totrep_frac]==expfrac_tot(i))).D2];
    H3=scatter(expfrac_tot(i),sum(points_late_D<points_late_D2)/length(points_late_D),'r','filled'); 
    axis([0 Inf 0 1])
    legend([H1,H2,H3],'Low rep. (0-25%)','Middle rep. (25-50%)','High rep. (50-75%)','Location','SouthWest');
    ylabel('Fraction fast process','FontSize',14,'FontName','Arial');
    xlabel('Replicated fraction of time point','FontSize',14,'FontName','Arial');
    title('Fraction of fibers following fast process by replicated fraction of time point');
end
im_paper1([sample_path '/FractFiber_FastProcess_TotRep_byRepFiber' sample_path],5,3);

%Each plot correspond to a time point;
figure;

for i=1:length(expfrac_tot)
    subplot(1,length(expfrac_tot),i)
 hold on
axis([0 Inf 0 Inf])
title(['Bin' num2str(i) ' (tot_f=' num2str(round(expfrac_tot(i)*100)/100) ')'],'FontSize',14,'FontName','Arial');
scatter([sample([sample.totrep_frac]==expfrac_tot(i)).D],[sample([sample.totrep_frac]==expfrac_tot(i)).D2],15,color_points,'filled'); 


%scatter(D{j}(D{j}<D2{j}),D2{j}(D{j}<D2{j}),'g','filled'); 
%scatter(D{j}(D{j}>D2{j}),D2{j}(D{j}>D2{j}),'r','filled');
% scatter(D{j}(T{j}==3),D2{j}(T{j}==3),'g','filled');
% scatter(D{j}(T{j}==4),D2{j}(T{j}==4),'k','filled');
% scatter(D{j}(T{j}==5),D2{j}(T{j}==5),'y','filled');
ylabel('Distance slow process','FontSize',14,'FontName','Arial');
xlabel('Distance fast process','FontSize',14,'FontName','Arial');
% if j==1
% legend('Cluster 1','Cluster 2');
% end
plot([0 1],[0 1],'color',color_line,'LineWidth',1.5);
hold off
end
im_paper1([sample_path '/Transition_TotRep_' sample_path],20,3);

%Fibers in bins divided by replicated fraction of the fiber;
%colors by time point
lower_bin_repfiber=[0    0.1071    0.2143    0.3214    0.4286    0.5357    0.6429];
upper_bin_repfiber=[0.1071    0.2143    0.3214    0.4286    0.5357    0.6429    0.7500];
figure;
for i=1:length(lower_bin_repfiber)
    subplot(1,7,i)
 hold on
axis([0 Inf 0 Inf])
title(['Bin' num2str(i) ' (fiber_f=' num2str(round((lower_bin_repfiber(i)+upper_bin_repfiber(i))*100)/200) ')'],'FontSize',14,'FontName','Arial');
points_early_D=[sample(([sample.totrep_frac]>0) & ([sample.totrep_frac]<=0.3) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D];
points_early_D2=[sample(([sample.totrep_frac]>0) & ([sample.totrep_frac]<=0.3) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D2];
H1=scatter(points_early_D,points_early_D2,'g','filled'); 
points_middle_D=[sample(([sample.totrep_frac]>0.3) & ([sample.totrep_frac]<=0.6) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D];
points_middle_D2=[sample(([sample.totrep_frac]>0.3) & ([sample.totrep_frac]<=0.6) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D2];
H2=scatter(points_middle_D,points_middle_D2,'y','filled'); 
points_late_D=[sample(([sample.totrep_frac]>0.6) & ([sample.totrep_frac]<=1) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D];
points_late_D2=[sample(([sample.totrep_frac]>0.6) & ([sample.totrep_frac]<=1) & ([sample.fracfiber]>lower_bin_repfiber(i)) & ([sample.fracfiber]<=upper_bin_repfiber(i))).D2];
H3=scatter(points_late_D,points_late_D2,'r','filled'); 

%scatter(D{j}(D{j}<D2{j}),D2{j}(D{j}<D2{j}),'g','filled'); 
%scatter(D{j}(D{j}>D2{j}),D2{j}(D{j}>D2{j}),'r','filled');
% scatter(D{j}(T{j}==3),D2{j}(T{j}==3),'g','filled');
% scatter(D{j}(T{j}==4),D2{j}(T{j}==4),'k','filled');
% scatter(D{j}(T{j}==5),D2{j}(T{j}==5),'y','filled');
ylabel('Distance slow process','FontSize',14,'FontName','Arial');
xlabel('Distance fast process','FontSize',14,'FontName','Arial');
if i==1
    legend('Early (0-30%)','Middle (30-60%)','Late (60-100%)');
end
plot([0 1],[0 1],'color',color_line,'LineWidth',1.5);
hold off
end
im_paper1([sample_path '/Transition_RepFiber_bytimepoint_' sample_path],20,3);

figure(9);
im_paper1([sample_path '/Transition_RepFiber_' sample_path],20,3);

figure(10);
im_paper1([sample_path '/Transition_RepFiber_' sample_path 'noclass'],20,3);

figure(5);
im_paper1([sample_path '/Single_CorrProfiles_' sample_path],20,3);

figure(6);
im_paper1([sample_path '/Average_CorrProfiles_' sample_path],24,3.5);

figure(11);
im_paper1([sample_path '/FractFiber_FastProcess_RepFiber_' sample_path],5,3);

close all