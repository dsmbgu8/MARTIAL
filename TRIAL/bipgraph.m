function [G,Gall] = bipgraph(A,B,thresh,verbose)
% Input:
%  A = [N x n] matrix of N pivot samples, each of dimensionality n
%  B = [N x n] matrix of N pivot samples, each of dimensionality n
%      each pair (A(i,:),B(i,:)) assumed to be in alignment (e.g., same
%      or similar classes)
%  thresh = excludes all i,j for which d(A(i,:),B(j,:)*T) > thresh
%  verbose = enable verbose output (set >1 for *very* verbose output)
% Output: 
%  G = sparse graph of distances between A and B less than thresh
%  Gall = graph of all distances betweeen A and B 

if nargin < 4
  verbose=0;
end

dmtx=sqrt(fast_dmtx(A',B'));

if verbose
  fprintf('bipgraph: min(dmtx): %f max(dmtx): %f thresh: %f (%d points below thresh)\n',...
          min(dmtx(:)),max(dmtx(:)),thresh,sum(dmtx(:)<=thresh))
end

if nargout == 2
  Gall = dmtx;
end

dmtx(dmtx>=thresh) = Inf;

G = sparse(dmtx);
if sum(G(:)==0)==numel(G)
  fprintf('bipgraph: no points below thresh %f\n',thresh)
  pause
end

