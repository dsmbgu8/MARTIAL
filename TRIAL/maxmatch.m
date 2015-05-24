function matchidx = maxmatch(G)
% Computes maximal matching for [N x M] bipartite graph G using hungarian algorithm
% Input:
%  G = N x M graph of all potential pairwise matches
% Output:
%  matchidx = N x 1 vector of matched indices 
%             (0 valued elements = not matched)
[matchidx,cost] = munkres(G); 
