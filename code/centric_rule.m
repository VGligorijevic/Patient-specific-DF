function Rnew = centric_rule(R,rule)
% Function for exporting significant values from the reconstructed matrix 
% by using row-centric or colomn-centric rule or combination (mix) of 
% these two
% -------------------------------------------------------------------------
% Vladimir Gligorijevic
% Imperial College London
% v.gligorijevic@imperial.ac.uk
% Last updated: 2/12/2014
% -------------------------------------------------------------------------
% [Input]: 
%   R: <2D Cell array>, r(node types) x r(node types) blocks
%   rule: <string>, type of export (e.g., row, col or mix})
% [Output]: 
%   R_new: <matrix>, matrix of significant values
%--------------------------------------------------------------------------

s = size(R);
if strcmp(rule,'row')
    row = repmat(mean(R,2),1,s(2));
    Rnew = R.*(R > row);
elseif strcmp(rule,'col')
    col = repmat(mean(R,1),s(1),1);
    Rnew = R.*(R  > col);
elseif strcmp(rule,'mix')
    col = repmat(mean(R,1),s(1),1);
    row = repmat(mean(R,2),1,s(2));
    Rnew = R.*((R > row).*(R  > col) > 0);
else
    fprintf('--Wrong type of export. Possible: {row, col, mix}\n');
end;
