function conn = connectivity(H)
% Function for creating connectivity matrix from matrix factor H
% -------------------------------------------------------------------------
% Vladimir Gligorijevic
% Imperial College London
% v.gligorijevic@imperial.ac.uk
% Last updated: 2/12/2014
% -------------------------------------------------------------------------
% [Input]:  
%   H: Cluster indicator matrix factor 
% [Output]: 
%   C: Binary connectivity matrix (Cij = 1 if pair (i,j) belong to same 
%      cluster, Cij = 0 otherwise)
% -------------------------------------------------------------------------

% If matrix is not properly transpose (first dimension must be larger, n >> k)
[n,k] = size(H);

% If matrix is not properly transpose (first dimension should be larger)
if k > n
    H = H';
end;

% Initialize connectivity matrix
conn = sparse(zeros(n,n));
[y,index]=max(H,[],2);  %find largest factor
mat1 = repmat(index',n,1);  % spread index down
mat2 = repmat(index,1,n); % spread index right

conn = mat1 == mat2; %connectivity matrix
