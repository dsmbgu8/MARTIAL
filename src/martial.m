function outparm = martial(S,SL,T,TL,SU,TU,Qk,options)
% Input: 
%  S = [NS x n] matrix of source samples, each of dimensionality n
%  SL = [NS x 1] vector of labels for source samples, each in [1,K]
%  T = [NT x n] matrix of labeled target samples, each of dimensionality n
%  TL = [NT x 1] vector of target labels (optional, use [] if unavailable)
%  SU = [MS x n] matrix of unlabeled source samples, each of dimensionality n
%  TU = [MT x n] matrix of unlabeled target samples, each of dimensionality n
%  Qk = number of pivots for each class 
%  options = struct of options, including:
%   options.Nrand = number of random initializations for seed
%             transformation (default 100)
%   options.Npool = number of pool pivot samples per class (default 300)
%   options.verbose = enable verbose output (set >1 for *very* verbose output)
% Output:
%  outparm = struct of output parameters, including:
%    outparm.SX[step] = source samples after applying 
%                       TRIAL (init|align|improve) step
%    outparm.PSXrmsd[step] = RMSD between pivots after applying
%                            TRIAL (init|align|improve) step 
%    outparm.PSXn[step] = number of pivots used in each of the 
%                         TRIAL (init|align|improve) step 
%    outparm.pivot_pool = struct of (at most) Npool*K pivot pool samples
%                         (see select_pivots.m for details)
%  
% NOTE: S, T, SU and TU assumed L^2 normalized 
%  (e.g., x=x/\|x\|_2 \forall x \in {S,T,SU,TU})

if ~exist('options','var')
  options = struct();
end

verbose=0;
if isfield(options,'verbose')
  verbose=options.verbose;
end

Nrand=100;
if isfield(options,'Nrand')
  Nrand=options.Nrand;
end

Npool=300;
if isfield(options,'Npool')
  Npool=options.Npool;
end

USL = unique(SL);
K = numel(USL);


% transformed samples for each TRIAL step
SXinit = [];
SXalign = [];
SXimprov = [];

% per class rmsd for each TRIAL step
PSXrmsdinit = [];
PSXrmsdalign = [];
PSXrmsdimprov = [];

% per class number of pivots for each TRIAL step
PSXninit = [];
PSXnalign = [];
PSXnimprov = [];


% pick a big pool of pivots
pivot_pool = select_pivots(S,SL,T,TL,Npool,SU,TU);
PSI=pivot_pool.srcTrueIdx;
PSL=pivot_pool.srcTrueLab;

PSTI=pivot_pool.tgt2srcIdx;
PSTL=pivot_pool.tgt2srcLab;

PSU=pivot_pool.srcUnlabIdx;
PSUL=pivot_pool.srcUnlabLab;

PSTU=pivot_pool.tgt2srcUnlabIdx;
PSTUL=pivot_pool.tgt2srcUnlabLab;


% merge the labeled and unlabeled pivots for filtering
PS = [S(PSI,:); SU(PSU,:)];
PSL = [PSL; PSUL];

PT = [T(PSTI,:); TU(PSTU,:)];
PTL = [PSTL; PSTUL];

% select best Npool source-to-target (labeled and unlabeled) pivots
do_filter=1; 
if do_filter
  [PSIPR,PTIPR] = filter_pivots(PS,PSL,PT,PTL,Npool);
  PS = PS(PSIPR,:);
  PSL = PSL(PSIPR);
  PT = PT(PTIPR,:);
  PTL = PTL(PTIPR);
end


for j = 1:K
  % get source samples for class j
  Sj = S(SL==USL(j),:);
  
  % get source and target pivots for class j
  PSidxj = find(PSL==USL(j));
  PTidxj = find(PTL==USL(j));
  PSj = PS(PSidxj,:);
  PTj = PT(PTidxj,:);
  NPj = numel(PSidxj);
  
  if NPj ~= numel(PTidxj)
    fprintf('martial: # source pivots (%d) != # tgt pivots (%d)',NPj,numel(PTidxj));
    pause
  end
  
  % compute initial transform using first Qk pivots
  PSinitj = PSidxj(1:min(Qk,NPj));
  PTinitj = PTidxj(1:min(Qk,NPj));
  [PSXformj,PSXtransj,PSXscalej,PSXrmsdj] = Kabsch(PT(PTinitj,:)',PS(PSinitj,:)');
  PSTrmsdbestj = PSXrmsdj; 
  if verbose>1
    dists = diag(fast_dmtx(PT(PTinitj,:)',PS(PSinitj,:)'));
    fprintf('martial: PSTrmsdbestj(%d) %0.4f, dists (%0.3f)=\n',j,PSTrmsdbestj,sum(dists)); fprintf('%0.3f ',dists); fprintf('\n')
    fprintf('martial: PSinitj= '); fprintf('%d ',PSinitj); fprintf('\n');
    fprintf('martial: PTinitj= '); fprintf('%d ',PTinitj); fprintf('\n');
  end

  randj = 0;
  while randj < Nrand
    % randomly select pivot pairs to find a better initial transform
    PSinitrj = PSidxj(randperm(NPj,Qk));
    PTinitrj = PTidxj(randperm(NPj,Qk));
        
    [PSXformrj,PSXtransrj,PSXscalerj,PSTrmsdrj] = Kabsch(PT(PTinitrj,:)',PS(PSinitrj,:)');            
    if PSTrmsdrj < PSTrmsdbestj
      PSinitj = PSinitrj;
      PTinitj = PTinitrj;
      PSXformj = PSXformrj;
      PSXtransj = PSXtransrj;
      PSTrmsdbestj = PSTrmsdrj;
      if verbose>1
        dists = diag(fast_dmtx(PT(PTinitj,:)',PS(PSinitj,:)'));
        fprintf('martial: PSTrmsdbestj(%d) %0.4f, dists (%0.3f)=\n',j,PSTrmsdbestj,sum(dists)); fprintf('%0.3f ',dists); fprintf('\n')
        fprintf('martial: PSinitj= '); fprintf('%d ',PSinitj); fprintf('\n');
        fprintf('martial: PTinitj= '); fprintf('%d ',PTinitj); fprintf('\n');
      end
    end
    randj = randj+1;      
  end
  [PSXforminitj,PSXtransj,PSXscalej,PSXrmsdj] = Kabsch(PT(PTinitj,:)',PS(PSinitj,:)');

  SXj = (Sj*PSXforminitj); %+repmat(PSXtransj, size(Sj,1), 1);
  SXinit = [SXinit; SXj];
  
  PSXrmsdinit = [PSXrmsdinit, PSXrmsdj];
  PSXninit = [PSXninit, numel(PTinitj)];
  
  %Smeanj = mean(Sj); % mean(PS(PSidxj,:));       
  threshj = PSXrmsdj; %rmsd(repmat(Smeanj, size(Sj,1), 1), Sj);  
  
  % update indices of PSinitj \in PTidxj and PTinitj \in PTidxj for align
  PSidxj = find(ismember(PSidxj,PSinitj));
  PTidxj = find(ismember(PTidxj,PTinitj));
  [PTidxj,PSidxj,PSXformalignj,PSXrmsdj] = align(PTidxj,PSidxj,PTj,PSj,...
                                               PSXformj,threshj,verbose);        

  SXj = (Sj*PSXformalignj); %+repmat(PSXtransj, size(Sj,1), 1);
  SXalign = [SXalign; SXj];
  PSXrmsdalign = [PSXrmsdalign, PSXrmsdj];
  PSXnalign = [PSXnalign, numel(PTidxj)];

  threshj = PSXrmsdj;
  [PTidxj,PSidxj,PSXformimprovj,PSXrmsdj] = improve(PTidxj,PSidxj,PTj,PSj,...
                                            PSXformalignj,threshj,verbose);        

  %SXj = (PSscalefj*(Sj-meanrepj)*PSXformj)+meanrepj;
  SXj = (Sj*PSXformimprovj); %+repmat(PSXtransj, size(Sj,1), 1);
  SXimprov = [SXimprov; SXj];
  PSXrmsdimprov = [PSXrmsdimprov, PSXrmsdj];
  PSXnimprov = [PSXnimprov, numel(PTidxj)];
end


outparm = struct();

outparm.SXinit = SXinit;
outparm.SXalign = SXalign;
outparm.SXimprov = SXimprov;

outparm.PSXrmsdinit = PSXrmsdinit;
outparm.PSXrmsdalign = PSXrmsdalign;
outparm.PSXrmsdimprov = PSXrmsdimprov;

outparm.PSXninit = PSXninit;
outparm.PSXnalign = PSXnalign;
outparm.PSXnimprov = PSXnimprov;

outparm.pivot_pool = pivot_pool;