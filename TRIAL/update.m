function [P,Q] = update(G,matchidx)
% G = [N x M] matrix representation of graph G 
% matchidx = [N x 1] vector of column indices (each in [1,M]) of matches
%            (0 valued elements = not matched)
% P = [Nmatched x 1] vector of row indices of matched samples
% Q = [Nmatched x 1] vector of column indices of matched samples 
matched = (matchidx~=0);
Nmatched = sum(matched);
P = 1:size(G,1);
P = P(matched);
Q = matchidx(matched);

