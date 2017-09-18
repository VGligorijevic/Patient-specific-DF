k = [3, 25, 20];
max_iter = 100;
initialization = 'random';

ext = [num2str(k(1)) '_' num2str(k(2)) '_' num2str(k(3))];

adj_list = {'./data/OVARIAN_PATIENT-PATIENT_mtrx.txt','./data/OVARIAN_GENE-GENE_mtrx.txt','./data/OVARIAN_DRUG-DRUG_mtrx.txt'};
rel_file = './data/OVARIAN_RELATIONS_mtrx.txt';
% create block matrices
[R, A, label_list] = block_matrices(adj_list, rel_file);

% run NMTF and export clusters
run_simNMTF(k, R, A, label_list, max_iter, initialization);

% compute survival times for patients in clusters
unix(['python survival.py ./data/OVARIAN_clinical.txt ./results/patients_final_' ext '_clust.txt  ./results/patients_final_' ext '.surv']);

% compute KM curves
D = kaplan_meier(['./results/patients_final_' ext '.surv']);