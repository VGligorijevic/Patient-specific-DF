function [S,G] = factorization_ssnmtf(R,A,k,max_iter,initialization)
% Function for simultaneous decomposition of all connections 
% -------------------------------------------------------------------------------------------------------------
% Vladimir Gligorijevic
% Imperial College London
% v.gligorijevic@imperial.ac.uk
% Last updated: 2/07/2015
% --------------------------------------------------------------------------------------------------------------
% Implementaiton based on the paper:
%
% Wang H., et al., Simultaneous  clustering  of  multi-type  relational  data  via  symmetric
% nonnegative matrix tri-factorization. In The 20th ACM Conference on
% Information and Knowledge Management (2011)
%
% [Input]:
%     R: <2D Cell array >, r(node types) x r(node types) blocks, relational matrix (e.g., R{i,j} = Rij, 
%     ni(nodes of type i) x nj(nodes of type j))
%     A: <1D Cell array>, r (node types) blocks, adjacency matrix (e.g., A{i} = Ai, ni (nodes) x ni (nodes)) 
%     k: <array>, rank parameters (e.g., [k1, k2,...,kr])
%     max_iter: <int>, predefined number of iterations 
%     initialization: <string>, initialization strategy: random, random_acol, nnmf
% [Output]: 
%     S: <2D Cell>, r(node types) x r(node types) blocks, compressed matrix (e.g., S{i,j} = Sij, ki x kj)
%     G: <1D Cell>, r(node types) blocks, cluster indicator matrix (e.g., G{i} = Gi, ni x ki)
% --------------------------------------------------------------------------------------------------------------

r = length(A);

% Compliting relation matrix
n = [];
for ii=1:r
    R{ii,ii} = A{ii};
    n(ii) = length(A{ii}); % sizes 
end;

% For of each data source
fprintf('-Initializations of G matrices....\n');
for ii=1:r
    % G matrix initialization
    G{ii} = matrix_initialization(R,ii,n(ii),k(ii),initialization); 
    fprintf('Initialization of G[%d] matrix finished!\n',ii);
end;
fprintf('-Initialization of G matrices finished!\n\n');


fprintf('-Initializations of S matrices....\n');
for ii=1:r
    for jj=1:r
        if (nnz(R{ii,jj}) ~= 0)
            S{ii,jj} = sparse(rand(k(ii),k(jj)));
            fprintf('Initialization of S[%d,%d] matrix finished!\n',ii,jj);
        else
            S{ii,jj} = sparse(zeros(k(ii),k(jj)));
        end;
    end;
end;
fprintf('-Initialization of S matrices finished!\n\n');


% Converting to matrix representation
R = cell2mat(R);
if (~isequal(R,R'))
    error('Relation matrix not symmetric!');
end;
G = blkdiag(G{:});
S = cell2mat(S);

% Norm (R)
norm_R = norm(R,'fro')^2;

J_old = 0; %initialization of objective function
%Iterations 
fprintf('| Iteration | Cost | RSE | Rel_Var_Cost | \n');
for iter=1:max_iter
    
     
    GtG = G'*G;
    GS = G*S;
    
    % first
    RGS = R*GS;
    GSGtGS = GS*GtG*S;
    
    % second
    GtRG = G'*R*G;
    GtGSGtG = GtG*S*GtG;
    
    % Update rules
    G =  G.*(RGS./(GSGtGS + 1e-8)).^0.25;
    S = S.*(GtRG./(GtGSGtG + 1e-8)); 
    
    % Computing the relative square error (RSE) every 10th iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % (n-1)th iteration 
    if mod(iter+1,10) == 0    
        % objective (cost) function
        J_old = (norm( R - G*S*G','fro'))^2;
    end;
    
    % nth iteration 
    if mod(iter,10) == 0    
        % objective (cost) function
        J_new = (norm( R - G*S*G','fro'))^2;
    
        
        % Errors
        RSE = J_new/norm_R;
        rel_var = abs(J_new - J_old)/abs(J_old);
        
        % Writing output
        fprintf('%d %0.5e %0.5e %0.5e\n', iter, J_new, RSE, rel_var);
    
        if (RSE <= 1e-3 | rel_var <= 1e-6) % checking for convergence
            break;
        end;
    
    end;

end; 

% Converiting back to block matrices (cells)
S = mat2cell(S,k,k);
GG = mat2cell(G,n,k);
G = {};
for i=1:r
    G{i} = GG{i,i};
end;
