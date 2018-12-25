function [y, energy] = occlusionInferenceTRW(pB, pC, edgeind, eadj)
%
% [y, energy] = occlusionInferenceTRW(pB, pC, adj)
%
% Solves for y given boundary and continuity likelihoods and an adjacency
% list for each edge.  The MRF is constructed as a directed graph (that is,
% Eij(0,1)~=Eij(1,0)).  Each edgelet in the image forms two nodes in the
% graph (one for occlusions in each direction).  A link between these two
% nodes disallows them to both be 1.  Links between adjacent nodes in the
% directed graph encode the conditional term P(ej | ei).  Auxialliary
% junction nodes disallow edglets to begin from nowhere or end nowhere.
%
% Inputs:
%  pB(nedglets, 2): unary likelihoods
%  pC(npairs): p(ei = 1 | ej=1)
%  edgind(npairs, 2): edgelet indices
%  eadj{nedglets}: edge adjacency 
%

% Initial:
%   E*_i = [ -log[P(ei=0 | x)]       0 ]
%   E*_ji = [ 0                         0       
%            -log[P(ei=0|ej=1, x)]      -log[P(ei=1|ej=1, x)] ]
%
% Convert to below if not on a border junction:
%   K_i = min_j(E*_ji(2,2))
%   E_i = [E*_i(0)      K_i] 
%   E_ji = [0                       0       
%           E*_ji(2,1)-E*_i(0)      E*_ji(2,2)-K_i ]
%
% Junctions:
%   Infinite penalty for contour end and start, except at border junctions
%   Small penalty for T-junctions


VL = 100000; % very large

nnodes = size(pB, 1) + size(pB, 1)*2;

j3 = zeros(1, 8, 'single');
j3(5) = VL;
j4 = single(1, 16, 'single');
j4(9) = VL;


%% Unary potentials

Ei_0 = single(-log((1-pB)+1E-10));
unary = cell(nnodes, 1);
for k = 1:size(pB, 1)
    unary{k} = [Ei_0(k) 0];
end 


%% Pairwise potentials

Eji_0 = single(-log((1-pC)+1E-10));
Eji_1 = single(-log(pC+1E-10));

% conditional energies


for k2 = 1:numel(juncts)
    if numel(juncts{k2}==3) 
        unary{k+k2} = j3;
    elseif numel(juncts{k2}==4)
        unary{k+k2} = j4;
    else
        error(['invalid junction size: ' num2str([k2 numel(juncts)])]);
    end
end





[lab, energy] = ...
    trwMinimizeEnergy(nlabels, edges, unaryPot, pairPot, niter, stopeps, doNotCheck)