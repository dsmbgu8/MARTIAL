function thresh_max = maxthresh(A,B,P,Q,T,thresh,verbose)
% Input:
%  P,Q = indices into A,B of aligned samples P,Q, respectively
%  A = [N x n] matrix of N pivot samples, each of dimensionality n
%  B = [N x n] matrix of N pivot samples, each of dimensionality n
%      each pair (A(i,:),B(i,:)) assumed to be in alignment (e.g., same
%      or similar classes)
%  T = [n x n] transformation matrix mapping samples in B to A
%  thresh = excludes all i,j for which d(A(i,:),B(j,:)*T) > thresh
%  verbose = enable verbose output (set >1 for *very* verbose output)
% Output:
%  thresh_max = maximum euclidean distance that will not increas rmsd threshold

[N,n] = size(A);
[M,n] = size(B);
K = numel(P);
if numel(Q) ~= K
  fprintf('maxthr: |P| = %s != |Q| = %s\n', num2str(numel(P)), num2str(numel(Q)));
  pause
end
D = rmsd(A(P,:),B(Q,:)*T);

thresh_max = sqrt(thresh^2 + ((K*(thresh^2 - D^2)) / (min(M,N)-K)));

if verbose
  fprintf('maxthresh: M: %d, N: %d, K: %d, D: %g, thresh: %f, maxthresh: %f\n',...
          M,N,K,D,thresh,thresh_max);
end