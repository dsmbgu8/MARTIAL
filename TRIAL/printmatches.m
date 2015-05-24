function printmatches(Gall,matchidx,thresh)
% Prints the match indices and weights according to thresh 
% Input:
%  Gall = graph of all potential matches
%  matchidx = indices into Gall of proposed matches
%  thresh = matches <= thresh considered acceptable matches

fprintf('%d matches with thresh %f, weight range=[%f,%f]\n',...
        sum(matchidx~=0),thresh,min(Gall(:)),max(Gall(:)));
for i=1:size(Gall,1)
  fprintf('\ti=%d ',i)          
  mi = matchidx(i);
  if mi==0
    [minv,mini] = min(Gall(i,:));
    fprintf('unmatched, nearest j=%d (weight=%f)\n',mini,minv)
  else
    fprintf('matched j=%d (weight=%f)\n',mi,Gall(i,mi))    
  end
end