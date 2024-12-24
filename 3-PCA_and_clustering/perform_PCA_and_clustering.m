
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['../2-Plot_autocorrelation/' sample_path '/valuecorrall_' sample_path '.mat']);
load(['../2-Plot_autocorrelation/' sample_path '/rf_' sample_path '.mat']);

num_bins=length(valuecorrall);
n = ceil(num_bins/2);

for i=1:num_bins
    
    f = valuecorrall{i};
    length_fibers=(length(f)+1)/2;
    f_sc=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %I consider the correlation functions between 0 and 25 kb for the
    %analysis; as the function is symmetrical, I take the reion from 0 to
    %25kb corresponding to the range [length_fibers:length_fibers+25] in
    %the matrix
    
    for j=1:length(f(:,length_fibers))
        f_sc(j,:) = f(j,length_fibers:length_fibers+25);
    end
    
    %I calculate the distance matrix
    D = pdist(f_sc,'correlation');
    
    %I plot the distance matrix
    figure(1); 
    haxis(i)=subplot(2,n,i);
    hold on;
    axis([-Inf Inf -Inf Inf])
    title(['Bin' num2str(i) ' (f=' num2str(round(rf(i)*100)/100) ')'],'FontSize',12,'FontName','Arial');
    imagesc(squareform(D)); 
    colorbar('location','southoutside','FontSize',10,'FontName','Arial','box','off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %I check if the data are normal
    % test_normality{i}=normalitytest(f_sc(:,3)');
    
    %PCA on f_sc
    [pc,score,latent,tsquare] = princomp(f_sc);
    var_exp=latent./sum(latent);
    figure(2); subplot(2,n,i);[hPareto, axesPareto] =pareto(var_exp); 
    yticks = get(axesPareto(2),'YTick');
    RightYLabels = cellstr(get(axesPareto(2),'YTickLabel'));
    set(axesPareto(1),'YTickLabel',[]);
    set(axesPareto(1),'YTickLabel',RightYLabels)
    set(axesPareto(2),'YTickLabel',[]);
    xlabel('Principal Component','fontname','Arial','fontsize',12);
    title(['Bin' num2str(i) ' (f=' num2str(round(rf(i)*100)/100) ')'],'FontSize',10,'FontName','Arial');
    %disp(['first + second components = ',num2str(sum(var_exp(1:2)))]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Silhouette comparison
    figure(3); subplot(2,n,i); hold on;
    title(['Bin' num2str(i) ' (f=' num2str(round(rf(i)*100)/100) ')'],'FontSize',10,'FontName','Arial');
    for n_clust=1:5
        %Hierarchical clustering 'weighted'
        Z = linkage(D,'weighted');
        T{i} = cluster(Z,'maxclust',n_clust);
        s = silhouette(f_sc,T{i},'correlation');
        silho_clust(1,n_clust)=nanmean(s);
    end
    scatter(1:5,silho_clust,'b','filled')
    plot(1:5,silho_clust,'b-','linewidth',1)
    ylabel('Averaged silhouette values','fontname','Arial','fontsize',12);
    xlabel('Number of clusters','fontname','Arial','fontsize',12);
    set(gca,'XLim',[1,6])
    set(gca,'XTick',[1,2,3,4,5,6])
    set(gca,'XTickLabel',{'1';'2';'3';'4';'5';'6'})
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Hierarchical clustering 'weighted'
    Z = linkage(D,'weighted');
    T{i} = cluster(Z,'maxclust',2);
    
    %I plot the number of fibers per bin
    figure(4);
    subplot(2,n,i);
    hold on
     title(['Bin' num2str(i) ' (f=' num2str(round(rf(i)*100)/100) ')'],'FontSize',12,'FontName','Arial');
    [a,b]=hist(T{i},1:2);
    axis([0.3 2.5 0 1])
    xlabel('Clusters','FontSize',12,'FontName','Arial');
    ylabel('Fraction of fibers','FontSize',12,'FontName','Arial');
    if a(1) ~= a(2)
        label_cluster_lessvalues=find(a == min(a));
        label_cluster_morevalues=find(a == max(a));
    else
        label_cluster_lessvalues=1;
        label_cluster_morevalues=2;
    end
    bar(label_cluster_lessvalues,a(label_cluster_lessvalues)/sum(a),'r');
    bar(label_cluster_morevalues,a(label_cluster_morevalues)/sum(a),'b');
    set(gca,'XTick',[1 2])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %I plot the averaged profiles of C(r,f) in the two classes
    fiber1= find(T{i} == label_cluster_lessvalues);
    fiber2 = find(T{i} == label_cluster_morevalues);
    figure(5);subplot(2,n,i);hold on;
    axis([-Inf Inf -Inf Inf])
    title(['Bin' num2str(i) ' (f=' num2str(round(rf(i)*100)/100) ')'],'FontSize',12,'FontName','Arial');
    ylabel('C(f,r) / C(f,0)','FontSize',12,'FontName','Arial');
    xlabel('r (kb)','FontSize',12,'FontName','Arial');
    if size(f_sc(fiber1,:),1) ~= 1
        errorbar([1: size(f_sc(fiber1,:),2)],mean(f_sc(fiber1,:)),std(f_sc(fiber1,:))/sqrt(length(fiber1)),'r','linewidth',1); 
    else
        plot([1: size(f_sc(fiber1,:),2)],f_sc(fiber1,:),'r','linewidth',1); 
    end 
    
    if size(f_sc(fiber2,:),1) ~= 1
        errorbar([1: size(f_sc(fiber2,:),2)],mean(f_sc(fiber2,:)),std(f_sc(fiber2,:))/sqrt(length(fiber2)),'b','linewidth',1); 
    else
        plot([1: size(f_sc(fiber2,:),2)],f_sc(fiber2,:),'b','linewidth',1);
    end
    legend('Cluster 1','Cluster 2')

end

save([sample_path '/T_' sample_path '.mat'],'T');

figure(1); 
hold on
for i=5:7
pos = get( haxis(i), 'Position' );
pos(2)=pos(2)+0.1;
set( haxis(i), 'Position',pos );
end
im_paper1([sample_path '/Corr_matrix_' sample_path ],10,6);

figure(2)
im_paper1([sample_path '/PCAPareto_' sample_path ],10,6);

figure(3)
im_paper1([sample_path '/Silhouette_' sample_path ],10,6);

figure(4);
im_paper1([sample_path '/Cluster_histogram_' sample_path ],10,6);

figure(5);
im_paper1([sample_path '/Cluster_meanprofile_' sample_path ],10,6);

close all