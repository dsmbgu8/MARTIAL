function val = rmsd(P,Q)
% Calculates the root-mean-squared deviation for K aligned samples P and Q
% Input:
%  P = [K x n] matrix 
%  Q = [K x n] matrix 
% Output:
%  val = Root-Mean-Square Deviation of P and Q
K = size(P,1);
if K ~= size(Q,1)
  fprintf('rmsd: |P| = %s != |Q| = %s\n', num2str(numel(P)), num2str(numel(Q)));
  pause
end
val = sqrt(sum(diag(fast_dmtx(P',Q')))/K);