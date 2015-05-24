function [P,Q,T,PQTrmsd] = align(P,Q,A,B,T,thresh,verbose) 
% Alignment step of TRIAL algorithm 
% Input:
%  (P,Q) = indices into (A,B) of aligned samples (P,Q), respectively
%  (A,B) = [N x n] matrices of N pivot samples, each of dimensionality n
%      each pair (A(i,:),B(i,:)) assumed to be in alignment (i.e., same
%      or similar classes)
%  T = [n x n] transformation matrix mapping samples in B to A
%  thresh = excludes all (i,j) for which d(A(i,:),B(j,:)*T) > thresh
%  verbose = enable verbose output (set >1 for *very* verbose output)
% Output: 
%  (P,Q) = indices into (A,B) of refined alignment
%  T = [n x n] transformation matrix mapping samples in B to A using
%      refined alignment
%  PQTrmsd = RMSD of refined alignment

if ~exist('verbose','var'), verbose=0; end;

[NA,n]=size(A);
[NB,n]=size(B);

if NA~=NB
  disp('align: A and B not aligned')
  pause
end

% handle transposed P and Q
if size(P,1) ~= 1, P = P'; end;
if size(Q,1) ~= 1, Q = Q'; end;

PQTrmsd = thresh; 
PQTrmsd_prev = thresh; 

Mp = numel(P); % number of matched samples, start with seed size

iter=0;
while(1)
  M = Mp;

  if M==numel(A), break; end; % all pivots in A matched to B

  % get matches w.r.t. threshold
  [G,Gall] = bipgraph(A,B*T,thresh);
  matchidx = maxmatch(G);
  
  if verbose>1
    fprintf('align[%d]: ',iter)
    printmatches(Gall,matchidx,thresh); 
  end

  Mp = sum(matchidx~=0);
  if (Mp > M) 
    [P,Q] = update(G,matchidx);
    
    if (numel(unique(P))~=numel(P)) | (numel(unique(Q))~=numel(Q))
      fprintf('align[%d]: non-unique P or Q!',iter)
      pause
    end
    
    T = Kabsch(A(P,:)',B(Q,:)');    

    PQTrmsd = rmsd(A(P,:),B(Q,:)*T);
  end
    
  if verbose
    fprintf('\talign[%d], M=%d, Mp=%d, PQRMSd=%g, PQRMSd_prev=%g\n',...
            iter,M,Mp,PQTrmsd,PQTrmsd_prev);
    PQTrmsd_prev = PQTrmsd;
  end
  iter = iter+1;
  if (Mp <= M), break; end
end


