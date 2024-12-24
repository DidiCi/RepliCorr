function storeexperfiber(unit,pag,path,filename,sample_path)
%Program to store the data from excel in variables to use in Matlab
%I store the data in a structure

for p=1:length(filename)
sprintf('Treatment of file n. %i',p)

%With "piece" I indicate the single fiber in the current file
num_pieces=xlsread([path filename{p}],pag,'B4'); %Number of fibers in the current file
length_pieces=xlsread([path filename{p}],pag,'D:D'); %Lengths of fibers in the current file

%I don't take anymore the replication values from the excel file because I
%calculate all from the fibers
%exfracfiber=xlsread(filename{p},pag,'F2:F110'); %Replication for fiber, F:F take the colomn
%extotfrac_rep=xlsread(filename{p},pag,'B3'); %Total fraction of replication for fibers

%I take from the excel file the data indicating replicated ad not replicated part,
%i convert them in string of 0 and 1 like my genome and I save in
%the structure exDcut ('ex' for experimental, 'D' for DNA, 'cut' because
%divided in fibers).
%The data from the excel file are round to the kb, so the length of the
%fibers and the replication level after the convertion could be a bit different
exDcut=[]; %I inizialize the structure
[exDcut(1:num_pieces).fiber] = deal([]); %I initialize all the fibers

for i=1:num_pieces
    temp=xlsread([path filename{p}],pag+i,'A:B'); 
    First_temp=xlsread([path filename{p}],pag+i,'A1'); 
    %If temp(1,1) is NaN the fiber start with a replicated part
    if isnan(temp(1,1))
    temp2=round(temp'*(1000/unit)); %Divided by 1000 because the unit in excel is 1kb
    temp2(isnan(temp2))=[];
    %I construct the O/1 fiber
    for j=1:length(temp2)
    if mod(j,2)~=0
        exDcut(i).fiber=[exDcut(i).fiber;ones(temp2(j),1)];
    end
    if mod(j,2)==0
        exDcut(i).fiber=[exDcut(i).fiber;zeros(temp2(j),1)];
    end
    end
    %If temp(1,1) is not NaN and First_temp is not empty the fiber start
    %with a  unreplicated part or is completely unreplicated
    elseif ~isnan(temp(1,1)) && ~isempty(First_temp)
    temp2=round(temp'*(1000/unit));
    temp2(isnan(temp2))=[];
    %I construct the O/1 fiber
    for j=1:length(temp2)
    if mod(j,2)==0
         exDcut(i).fiber=[exDcut(i).fiber;ones(temp2(j),1)];
    end
    if mod(j,2)~=0
         exDcut(i).fiber=[exDcut(i).fiber;zeros(temp2(j),1)];
    end
    end
    elseif ~isnan(temp(1,1)) && isempty(First_temp)
        temp2=round(temp'*(1000/unit));
        temp2(isnan(temp2))=[];
        %I construct the O/1 fiber
        exDcut(i).fiber=[exDcut(i).fiber;ones(temp2,1)];
    end
    
end

%I indicate for each fiber the unit of the sigle block in bp
[exDcut.unit_block]=deal(unit);

%I create a bigger structure with all the fibers of all the files to
%analize, another structure for the number of fibers in each file and
%another for the length of the fibers in each file
globalallexDcut.(['exDcut' filename{p}])=exDcut;
globalallnum_pieces.(['num_pieces' filename{p}])=num_pieces
globalalllength_pieces.(['length_pieces' filename{p}])=length_pieces;


save([sample_path '/globalallexDcut.mat'],'globalallexDcut') 
save([sample_path '/globalallnum_pieces.mat'],'globalallnum_pieces')
save([sample_path '/globalalllength_pieces.mat'],'globalalllength_pieces') 
end


end