function [lab, energy, labs, tmpim] = occlusionInferenceLBP(pB, pC, edgeind, bndinfo)
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



%% Initialize edge data

ne = bndinfo.ne;

% unary
%pB = [pB(:, 1)  pB(:, 2) ; ...
%      pB(:, 1)  pB(:, 3)];
tic

% minimum energy (or maximum probability) for each directed edgelet being 
% an occluding contour
% put this energy into unary and divide it out of other terms
maxPon = pB(:, 2);
for k = 1:size(edgeind, 1)
    maxPon(edgeind(k, 2)) = max(maxPon(edgeind(k, 2)), pC(k));
end
for k = 1:size(edgeind, 1)
    pC(k) = pC(k) ./ maxPon(edgeind(k, 2));    
end

% continuity likelihood matrix
pCm = spalloc(ne*2, ne*2, size(edgeind, 1));
pCm(edgeind(:,2) + ne*2*(edgeind(:, 1)-1)) = pC; % P(e_out=1 | e_in=1, x)
% for k = 1:size(edgeind, 1)
%     pCm(edgeind(k, 2), edgeind(k, 1)) = pC(k);  % P(e_out=1 | e_in=1, x)
% end
toc

% relative angles
ra = spalloc(ne*2, ne*2, size(edgeind, 2));
theta = bndinfo.edges.thetaDirected*180/pi;
theta = [theta ; theta + 180];
theta = mod(theta, 360);
eadj = bndinfo.edges.adjacency;
for k = 1:numel(eadj)  
    for k2 = eadj{k}
        ra(k, k2) = mod(theta(k2)-theta(k), 360);
    end
end

%% Initialize junction data

nj = bndinfo.nj;
jund = repmat({[0 0 0 0]}, nj, 1);
jin = repmat({[0 0 0 0]}, nj, 1);
jout = repmat({[0 0 0 0]}, nj, 1);
jsize = zeros(nj, 1);
ejunctions = bndinfo.edges.junctions;
for k = 1:size(ejunctions, 1)
    j1 = ejunctions(k, 1);
    j2 = ejunctions(k, 2);
        
    jsize(j1) = jsize(j1)+1;
    jsize(j2) = jsize(j2)+1;
    
    jout{j1}(jsize(j1)) = k;  % e_k: j1 --> j2
    jin{j2}(jsize(j2)) = k;
    jout{j2}(jsize(j2)) = k+ne; % e_(k+ne): j2-->j1
    jin{j1}(jsize(j1)) = k+ne;
    jund{j1}(jsize(j1)) = k;
    jund{j2}(jsize(j2)) = k;
end
for k = 1:nj
    [jund{k}, jind] = sort(jund{k}(1:jsize(k)));
    jin{k} = jin{k}(1:jsize(k));
    jout{k} = jout{k}(1:jsize(k));
    jin{k} = jin{k}(jind);
    jout{k} = jout{k}(jind);
end


%% Initialize factor graph
nfactors = ne + nj;
factorLookup = spalloc(ne, nfactors, ne+sum(jsize));
factors = cell(nfactors, 1);
nodeSizes = repmat(3, ne, 1);

%% Unary factors
sz = 3;
for k = 1:ne
    factorLookup(k, k) = 1;
    jpot = [pB(k, 1) ;  maxPon(k) ; maxPon(ne+k)];
    factors{k} = tabular_kernel(sz, jpot);
end 

tic
%% Junction factors
for k = 1:nj
    
    fk = k + ne;    
    
    switch jsize(k)
        
        case 1            
            out1 = 2+(jout{k}(1)>ne); in1 = (4-out1)+1;  % out and in = 2 or 3
            jpot = 1E-10*ones(3,1);
            sz = 3;            
            
            jpot(1) = 1;
            jpot(out1) = pB(jout{k}, 2);
            jpot(in1) = 1;                                   
            
        case 2      
            out1 = 2+(jout{k}(1)>ne); in1 = (4-out1)+1;  % out and in = 2 or 3
            out2 = 2+(jout{k}(2)>ne); in2 = (4-out2)+1;
            jpot = 1E-10*ones(3,3);
            sz = [3 3];
            
            jpot(1, 1) = 1;  
            jpot(out1, in2) = pCm(jout{k}(1), jin{k}(2));
            jpot(in1, out2) = pCm(jout{k}(2), jin{k}(1));
        
        case 3
            ji1 = jin{k}(1);  ji2 = jin{k}(2);  ji3 = jin{k}(3);
            jo1 = jout{k}(1);  jo2 = jout{k}(2);  jo3 = jout{k}(3);
            out1 = 2+(jo1>ne); in1 = (4-out1)+1;  
            out2 = 2+(jo2>ne); in2 = (4-out2)+1;
            out3 = 2+(jo3>ne); in3 = (4-out3)+1;  
            jpot = 1E-10*ones(3,3,3);
            sz = [3 3 3];
            
            % zero in, zero out
            jpot(1,1,1) = 1;
            
            % one in, one out
            jpot(out1, in2, 1) = pCm(jo1, ji2);
            jpot(out1, 1, in3) = pCm(jo1, ji3);
            jpot(in1, out2, 1) = pCm(jo2, ji1);
            jpot(1, out2, in3) = pCm(jo2, ji3);
            jpot(in1, 1, out3) = pCm(jo3, ji1);
            jpot(1, in2, out3) = pCm(jo3, ji2);                                       
            
            % two in, one out: ra = rel angle (larger means more leftward) 
            jpot(out1, in2, in3) = ...
                (ra(ji2, jo1)>ra(ji3, jo1))*pCm(jo1, ji2) + ...
                (ra(ji2, jo1)<=ra(ji3, jo1))*pCm(jo1, ji3);
            jpot(in1, out2, in3) = ...
                (ra(ji1, jo2)>ra(ji3, jo2))*pCm(jo2, ji1) + ...
                (ra(ji1, jo2)<=ra(ji3, jo2))*pCm(jo2, ji3);            
            jpot(in1, in2, out3) = ...
                (ra(ji1, jo3)>ra(ji2, jo3))*pCm(jo3, ji1) + ...
                (ra(ji1, jo3)<=ra(ji2, jo3))*pCm(jo3, ji2);            
            
            % one in, two out 
            jpot(out1, out2, in3) = ...
                (ra(ji3, jo1)>ra(ji3, jo2))*pCm(jo1, ji3)*pB(jo2, 2) + ...
                (ra(ji3, jo1)<=ra(ji3, jo2))*pCm(jo2, ji3)*pB(jo1, 2);
            jpot(out1, in2, out3) = ...
                (ra(ji2, jo1)>ra(ji2, jo3))*pCm(jo1, ji2)*pB(jo3, 2) + ...
                (ra(ji2, jo1)<=ra(ji2, jo3))*pCm(jo3, ji2)*pB(jo1, 2);            
            jpot(in1, out2, out3) = ...
                (ra(ji1, jo2)>ra(ji1, jo3))*pCm(jo2, ji1)*pB(jo3, 2) + ...
                (ra(ji1, jo2)<=ra(ji1, jo3))*pCm(jo3, ji1)*pB(jo2, 2);   
            
        case 4
            ji1 = jin{k}(1);  ji2 = jin{k}(2);  ji3 = jin{k}(3);  ji4 = jin{k}(4);
            jo1 = jout{k}(1);  jo2 = jout{k}(2);  jo3 = jout{k}(3);  jo4 = jout{k}(4);
            out1 = 2+(jo1>ne); in1 = (4-out1)+1;  
            out2 = 2+(jo2>ne); in2 = (4-out2)+1;
            out3 = 2+(jo3>ne); in3 = (4-out3)+1;  
            out4 = 2+(jo4>ne); in4 = (4-out4)+1; 
            jpot = 1E-10*ones(3,3,3,3);
            sz = [3 3 3 3];
            
            % zero in, zero out
            jpot(1,1,1,1) = 1;
            
            % one in, one out
            jpot(out1, in2, 1, 1) = pCm(jo1, ji2);           
            jpot(out1, 1, in3, 1) = pCm(jo1, ji3);
            jpot(out1, 1, 1, in4) = pCm(jo1, ji4);            
            jpot(in1, out2, 1, 1) = pCm(jo2, ji1);
            jpot(1, out2, in3, 1) = pCm(jo2, ji3);
            jpot(1, out2, 1, in4) = pCm(jo2, ji4);
            jpot(in1, 1, out3, 1) = pCm(jo3, ji1);
            jpot(1, in2, out3, 1) = pCm(jo3, ji2);                                       
            jpot(1, 1, out3, in4) = pCm(jo3, ji4);
            jpot(in1, 1, 1, out4) = pCm(jo4, ji1);
            jpot(1, in2, 1, out4) = pCm(jo4, ji2);                                       
            jpot(1, 1, in3, out4) = pCm(jo4, ji3);
            
            % two in, one out: ra = rel angle (larger means more leftward) 
            jpot(out1, in2, in3, 1) = ...
                (ra(ji2, jo1)>ra(ji3, jo1))*pCm(jo1, ji2) + ...
                (ra(ji2, jo1)<=ra(ji3, jo1))*pCm(jo1, ji3);
            jpot(out1, in2, 1, in4) = ...
                (ra(ji2, jo1)>ra(ji4, jo1))*pCm(jo1, ji2) + ...
                (ra(ji2, jo1)<=ra(ji4, jo1))*pCm(jo1, ji4);
            jpot(out1, 1, in3, in4) = ...
                (ra(ji3, jo1)>ra(ji4, jo1))*pCm(jo1, ji3) + ...
                (ra(ji3, jo1)<=ra(ji4, jo1))*pCm(jo1, ji4);            
            jpot(in1, out2, in3, 1) = ...
                (ra(ji1, jo2)>ra(ji3, jo2))*pCm(jo2, ji1) + ...
                (ra(ji1, jo2)<=ra(ji3, jo2))*pCm(jo2, ji3);            
            jpot(in1, out2, 1, in4) = ...
                (ra(ji1, jo2)>ra(ji4, jo2))*pCm(jo2, ji1) + ...
                (ra(ji1, jo2)<=ra(ji4, jo2))*pCm(jo2, ji4);  
            jpot(1, out2, in3, in4) = ...
                (ra(ji3, jo2)>ra(ji4, jo2))*pCm(jo2, ji3) + ...
                (ra(ji3, jo2)<=ra(ji4, jo2))*pCm(jo2, ji4);              
            jpot(in1, in2, out3, 1) = ...
                (ra(ji1, jo3)>ra(ji2, jo3))*pCm(jo3, ji1) + ...
                (ra(ji1, jo3)<=ra(ji2, jo3))*pCm(jo3, ji2); 
            jpot(in1, 1, out3, in4) = ...
                (ra(ji1, jo3)>ra(ji4, jo3))*pCm(jo3, ji1) + ...
                (ra(ji1, jo3)<=ra(ji4, jo3))*pCm(jo3, ji4); 
            jpot(1, in2, out3, in4) = ...
                (ra(ji2, jo3)>ra(ji4, jo3))*pCm(jo3, ji2) + ...
                (ra(ji2, jo3)<=ra(ji4, jo3))*pCm(jo3, ji4);                         
            jpot(in1, in2, 1, out4) = ...
                (ra(ji1, jo4)>ra(ji2, jo4))*pCm(jo4, ji1) + ...
                (ra(ji1, jo4)<=ra(ji2, jo4))*pCm(jo4, ji2); 
            jpot(in1, 1, in3, out4) = ...
                (ra(ji1, jo4)>ra(ji3, jo4))*pCm(jo4, ji1) + ...
                (ra(ji1, jo4)<=ra(ji3, jo4))*pCm(jo4, ji3); 
            jpot(1, in2, in3, out4) = ...
                (ra(ji2, jo4)>ra(ji3, jo4))*pCm(jo4, ji2) + ...
                (ra(ji2, jo4)<=ra(ji3, jo4))*pCm(jo4, ji3);
            
            % one in, two out 
            jpot(in1, out2, out3, 1) = ...
                (ra(ji1, jo2)>ra(ji1, jo3))*pCm(jo2, ji1)*pB(jo3, 2) + ...
                (ra(ji1, jo2)<=ra(ji1, jo3))*pCm(jo3, ji1)*pB(jo2, 2);
            jpot(in1, out2, 1, out4) = ...
                (ra(ji1, jo2)>ra(ji1, jo4))*pCm(jo2, ji1)*pB(jo4, 2) + ...
                (ra(ji1, jo2)<=ra(ji1, jo4))*pCm(jo4, ji1)*pB(jo2, 2);   
            jpot(in1, 1, out3, out4) = ...
                (ra(ji1, jo3)>ra(ji1, jo4))*pCm(jo3, ji1)*pB(jo4, 2) + ...
                (ra(ji1, jo3)<=ra(ji1, jo4))*pCm(jo4, ji1)*pB(jo3, 2);            
            jpot(out1, in2, out3, 1) = ...
                (ra(ji2, jo1)>ra(ji2, jo3))*pCm(jo1, ji2)*pB(jo3, 2) + ...
                (ra(ji2, jo1)<=ra(ji2, jo3))*pCm(jo3, ji2)*pB(jo1, 2); 
            jpot(out1, in2, 1, out4) = ...
                (ra(ji2, jo1)>ra(ji2, jo4))*pCm(jo1, ji2)*pB(jo4, 2) + ...
                (ra(ji2, jo1)<=ra(ji2, jo4))*pCm(jo4, ji2)*pB(jo1, 2);  
            jpot(1, in2, out3, out4) = ...
                (ra(ji2, jo3)>ra(ji2, jo4))*pCm(jo3, ji2)*pB(jo4, 2) + ...
                (ra(ji2, jo3)<=ra(ji2, jo4))*pCm(jo4, ji2)*pB(jo3, 2);                          
            jpot(out1, out2, in3, 1) = ...
                (ra(ji3, jo1)>ra(ji3, jo2))*pCm(jo1, ji3)*pB(jo2, 2) + ...
                (ra(ji3, jo1)<=ra(ji3, jo2))*pCm(jo2, ji3)*pB(jo1, 2);
            jpot(out1, 1, in3, out4) = ...
                (ra(ji3, jo1)>ra(ji3, jo4))*pCm(jo1, ji3)*pB(jo4, 2) + ...
                (ra(ji3, jo1)<=ra(ji3, jo4))*pCm(jo4, ji3)*pB(jo1, 2);
            jpot(1, out2, in3, out4) = ...
                (ra(ji3, jo2)>ra(ji3, jo4))*pCm(jo2, ji3)*pB(jo4, 2) + ...
                (ra(ji3, jo2)<=ra(ji3, jo4))*pCm(jo4, ji3)*pB(jo2, 2);                         
            jpot(out1, out2, 1, in4) = ...
                (ra(ji4, jo1)>ra(ji4, jo2))*pCm(jo1, ji4)*pB(jo2, 2) + ...
                (ra(ji4, jo1)<=ra(ji4, jo2))*pCm(jo2, ji4)*pB(jo1, 2);
            jpot(out1, 1, out3, in4) = ...
                (ra(ji4, jo1)>ra(ji4, jo3))*pCm(jo1, ji4)*pB(jo3, 2) + ...
                (ra(ji4, jo1)<=ra(ji4, jo3))*pCm(jo3, ji4)*pB(jo1, 2);   
            jpot(1, out2, out3, in4) = ...
                (ra(ji4, jo2)>ra(ji4, jo3))*pCm(jo2, ji4)*pB(jo3, 2) + ...
                (ra(ji4, jo2)<=ra(ji4, jo3))*pCm(jo3, ji4)*pB(jo2, 2);                                    
        otherwise
            error('junction greater than 4')
    end
    
    factorLookup(jund{k}, fk) = 1;
    factors{fk} = tabular_kernel(sz, jpot);

end
toc


%% Run inference

fg = mk_fgraph(factorLookup, nodeSizes, factors);

% get ground truth
gtlab = bndinfo.edges.boundaryType;
gtlab = 1 + (gtlab(1:end/2)>0) + (gtlab(end/2+1:end)>0)*2;
tmpim = zeros(size(bndinfo.wseg));
for k = find(gtlab>1)'   
    tmpim(bndinfo.edges.indices{k}) = 1;
end
figure(1), imagesc(tmpim), axis image, colormap gray
energy =  computeEnergy(fg.dom, factors, gtlab);
disp(['Ground Truth Energy: ' num2str(energy)]);

maximize = 0;
engine = belprop_fg_inf_engine(fg, 'maximize', maximize, 'max_iter', 1000);

evidence = cell(ne, 1);
[engine, loglikelihood, labs, niter] = enter_evidence_LowMem(engine, evidence, 'maximize', maximize);

lab = zeros(ne, 1);
for k = 1:ne
    m = marginal_nodes(engine, k);
    lab(k) = argmax(m.T);
end

energy = computeEnergy(fg.dom, factors, lab);

disp(['Error: ' num2str(mean(lab==gtlab))])
disp(['Energy: ' num2str(energy)]);

% convert to directed format
lab = [(lab==2) ; (lab==3)];

tmpim = zeros(size(bndinfo.wseg));
for k = 1:ne
    m = marginal_nodes(engine, k); 
    tmpim(bndinfo.edges.indices{k}) = 1-m.T(1);
end
figure(2), imagesc(tmpim), axis image
   

