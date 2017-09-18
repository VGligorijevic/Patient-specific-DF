% Test ranks
function run_simNMTF(k, R, A, label_list, max_iter, initialization)

ext = [num2str(k(1)) '_' num2str(k(2)) '_' num2str(k(3))];
[S, G] = factorization_ssnmtf(R,A,k,max_iter,initialization);
file_list = {['./results/patients_final_' ext],...
             ['./results/genes_final_' ext],...
             ['./results/drugs_final_' ext]};

% Exporting clusters
compute_clusters_ssnmtf(G, label_list, file_list);

% Writing connectivity matrix into file
dlmwrite(['./results/gclust-gclust_' ext '.mtrx'], full(S{2,2}), 'delimiter', ',');

% Exporting predictions
export_significant_associations(G, S, label_list, [2,2], ['./results/gene-gene_pred_' ext '.txt'], 'mix')
export_significant_associations(G, S, label_list, [2,3], ['./results/gene-drug_pred_' ext '.txt'], 'mix')
