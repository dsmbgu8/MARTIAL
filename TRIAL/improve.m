function [P,Q,T,PQTrmsd] = improve(P,Q,A,B,T,thresh,verbose)  
% Iterative improvement step of TRIAL algorithm 
% Input:
%  (P,Q) = indices into (A,B) of aligned samples (P,Q), respectively
%  (A,B) = [N x n] matrices of N pivot samples, each of dimensionality n
%      each pair (A(P(i),:),B(Q(i),:)) assumed to be in alignment (i.e., same
%      or similar classes)
%  T = [n x n] transformation matrix mapping samples in B to A
%  thresh = excludes all i,j for which d(A(i,:),B(j,:)*T) > thresh
%  verbose = enable verbose output (set >1 for *very* verbose output)
% Output: 
%  (P,Q) = indices into (A,B) of refined alignment
%  T = [n x n] transformation matrix mapping samples in B to A using
%      refined alignment
%  PQTrmsd = RMSD of refined alignment



if ~exist('verbose', 'var'), verbose=0; end;

Mp = numel(P); % initial set of matches
[NA,n]=size(A);
[NB,n]=size(B);

if Mp == NA
  if verbose, fprintf('improve[0,0]: all points already matched\n'); end;
  [T,PQTrmsd] = Kabsch(A(P,:)',B(Q,:)');
  return
end

% initial set of indices
IA = 1:NA;
IB = 1:NB;

% handle transposed P and Q
if size(P,1) ~= 1, P = P'; end;
if size(Q,1) ~= 1, Q = Q'; end;

iter0=0;
iter1=0;
while(1)
  Bp = B*T;
  M = Mp;
  
  PQTrmsd = rmsd(A(P,:),Bp(Q,:));
  if verbose>1
    fprintf('improve[%d,%d]: PQTrmsd before adding samples: %f (thresh=%f)\n',...
            iter0,iter1,PQTrmsd,thresh); 
  end
  
  while(1)
    thresh_max = maxthresh(A,B,P,Q,T,thresh,verbose);    
    if thresh_max == Inf | thresh_max <= eps
      if verbose
        fprintf('improve[%d,%d]: infinite or trivial threshold\n',iter0,iter1); 
      end
      break
    end
    
    % get samples \in A not currently \in P
    AminPI = setdiff(IA,P); 
    AminP = A(AminPI,:);
    
    % get samples \in B not currently \in Q
    BminQI = setdiff(IB,Q); 
    BminQ = B(BminQI,:);    

    [G,Gall] = bipgraph(AminP,BminQ,thresh_max);
    matchidx = maxmatch(G);
    if verbose>1
      fprintf('improve[%d,%d]: ',iter0,iter1)
      printmatches(Gall,matchidx,thresh); 
    end

    Mpp = sum(matchidx~=0);
    if (Mpp > 0)
      % update number and set of matched samples
      Mp = Mp+Mpp; 
      [Pm,Qm] = update(G,matchidx);

      P = [P AminPI(Pm)];
      Q = [Q BminQI(Qm)];
      
      if (numel(unique(P))~=numel(P)) | (numel(unique(Q))~=numel(Q))
        fprintf('improve[%d,%d]: non-unique P or Q!',iter0,iter1)
        pause
      end
      
    end
    PQTrmsd = rmsd(A(P,:),Bp(Q,:));
    if verbose
      fprintf('improve[%d,%d]: M=%d, Mp=%d, Mpp=%d, PQTrmsd(P,Q*T)=%g\n',...
              iter0,iter1,M,Mp,Mpp,PQTrmsd)
    end
    iter1 = iter1+1;
    if (Mpp <= 0), break; end
  end
  % recompute transform with larger P and Q
  T = Kabsch(A(P,:)',B(Q,:)');  
  PQTrmsd = rmsd(A(P,:),B(Q,:)*T);
  if verbose>1
    fprintf('improve[%d,%d]: PQTrmsd after adding samples: %f\n',...
            iter0,iter1,PQTrmsd); 
  end
  if (numel(P)~=numel(unique(P))) | (numel(Q)~=numel(unique(Q))) 
    fprintf('improve[%d,%d]: WARNING #P=%d, #uniq(P)=%d, #Q=%d, #uniq(Q)=%d\n',...
            iter0, iter1, numel(P),numel(unique(P)),numel(Q),numel(unique(Q)));
    pause
  end

  iter0 = iter0+1;
  if (Mp <= M), break; end;
end
