% Function for assigning entities to clusters
% -------------------------------------------------------------------------------------------------------------
% Vladimir Gligorijevic
% Imperial College London
% v.gligorijevic@imperial.ac.uk
% Last updated: 5/07/2015
% --------------------------------------------------------------------------------------------------------------
% [Input]:
%     G: <Cell list>, Cell list containing cluster indicator matrices
%     label_list: <Cell string>, arrays of labels for each data source
%     file_list: <Cell string>, list of filenames for exporting figures and
%     clusters
% --------------------------------------------------------------------------------------------------------------

function compute_clusters_ssnmtf(G, label_list, file_list)
fprintf('################################\n');
fprintf('Computing cluster assignment....\n');

for i=1:length(G)

    % Cluster indicator matrix
    [n,k] = size(G{i});

    if k > n
        G{i} = G{i}';
    end;

    [y,index] = max(G{i},[],2);  %find largest factor in column
    % Computing connectivity matrix
    C{i} = connectivity(G{i});


    % Exporting cluster indices into file
    fWrite = fopen([file_list{i} '_clust.txt'],'w');
    for ii=1:length(index)
        fprintf(fWrite,'%s %d\n',label_list{i}{ii},index(ii));
    end;
    fclose(fWrite);

    % Writing connectivity matrix into file
    dlmwrite([file_list{i} '.mtrx'],C{i},'delimiter',',');
    fprintf('Dataset [%d] finished.\n',i);
end;
