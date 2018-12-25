function alpha = trainEdgesPL2(bndinfo, pB, pC, edgeind, Y)
% alpha = trainEdgesPL(bndinfo, dtBnd, dtCont, doAll)
% P(e | x) = 1/Z * prod_i[ f0(e_i) f1(e_i, x) f2(e_i, parent(i), x) ]
% 
% doAll specifies whether parent(i) is only true for edges that are on in
% the ground truth

for f = 1:numel(bndinfo)

    %disp(num2str(f))
    
    ne = bndinfo(f).ne;
    nde = ne*2;
    
    % get unary term                                          
    x1{f} = pB{f};             
    
    isset2 = false(nde, 1);
    x2{f} = zeros(nde, 1);             
    
    spLR = bndinfo(f).edges.spLR;
    spLR = [spLR ; spLR(:, [2 1])];
    for k = 1:size(edgeind{f}, 1)
        k1 = edgeind{f}(k, 1);
        k2 = edgeind{f}(k, 2);             
        lab1 = bndinfo(f).labels(spLR(k1, 1));
        lab2 = bndinfo(f).labels(spLR(k2, 1));            
        if Y{f}(k1) > 0 && lab1==lab2
            x2{f}(k2) = pC{f}(k); % / tmppB(k2, 2); % / dtBnd.prior(1+(k2>ne));
            isset2(k2) = true;
        end        
    end  
%     x2{f} = [1-x2{f} x2{f}];
%     x2{f}(~isset2, :) = 1;        
    
    x2{f}(~isset2) = pB{f}(~isset2, 2);  
    x2{f} = [1-x2{f} x2{f}];
    
%    x2{f} = [ones(size(x2{f})) x2{f}];
%    x2{f}(~isset2, 2) = pB{f}(~isset2, 2);  

%     prior{f} = ones(nde, 2);
%     prior{f}(1:nde, 1) = dtBnd.prior(1);
%     prior{f}(1:ne, 2) = dtBnd.prior(2);
%     prior{f}(ne+1:nde, 2) = dtBnd.prior(3);
    
    ind = (Y{f}(1:ne)==-1) | (Y{f}(ne+1:end)==-1);
    ind = [ind ; ind];        
    
    ty{f} = Y{f}>0;
    
    x1{f}(ind, :) = [];
    x2{f}(ind, :) = [];
    %prior{f}(ind, :) = [];
    ty{f}(ind) = [];
    

end

%xp = vertcat(prior{:});
x1 = vertcat(x1{:});
x2 = vertcat(x2{:});
ty = vertcat(ty{:})';

ndata = numel(ty);

x{1} = [zeros(ndata, 1)  ones(ndata, 1)]; %xp; %repmat(dtBnd.prior, [ndata 1]);
x{2} = log(x1);
x{3} = log(x2);

%save './data/tmp.mat' x ty

% else
%     load './data/tmp.mat';
% end

alpha = [0 1 0];
options = optimset('Display', 'iter');
%alpha1 = fminsearch(@(p) negloglikelihood(p, x, ty), alpha, options)
alpha = fminunc(@(p) negloglikelihood(p, x, ty), alpha, options);
pval = getlikelihood(alpha, x, ty);
err = mean(pval < 0.5)
err0 = mean(pval(ty==0)<0.5)
err1 = mean(pval(ty==1)<0.5)

keyboard;


%% maximize likelihood (minimize -log likelihood)
%P(e_i=k|x, e) propto a(1)*x{1}(i,k) + a(2)*x{2}(i,k) + a(3)*x{3}(i,k)
function nll = negloglikelihood(alpha, x, y)
ndata = numel(y);
logpx = zeros(ndata, size(x{1}, 2));
for k = 1:numel(alpha)
    logpx = logpx + alpha(k)*x{k};
end
px = exp(logpx);
sumpx = sum(px, 2);
px = px ./ repmat(sumpx, [1 size(px, 2)]);

correctind = (1:ndata) + y*ndata;
pval = px(correctind);
disp(['ave p: ' num2str([sum(pval)/ndata sum(pval(y==0))/sum(y==0) ...
    sum(pval(y==1))/sum(y==1)])]); %  sum(pval(y==2))/sum(y==2)])]);
nll = -log(pval);
nll = sum(nll) / ndata;



%% get likelihood
function [pval, px] = getlikelihood(alpha, x, y)
ndata = numel(y);
logpx = zeros(ndata, size(x{1}, 2));
for k = 1:numel(alpha)
    logpx = logpx + alpha(k)*x{k};
end
px = exp(logpx);
sumpx = sum(px, 2);
px = px ./ repmat(sumpx, [1 size(px, 2)]);

correctind = (1:ndata) + y*ndata;
pval = px(correctind);
% disp(['ave p: ' num2str([sum(pval)/ndata sum(pval(y==0))/sum(y==0) ...
%     sum(pval(y==1))/sum(y==1)])]); %  sum(pval(y==2))/sum(y==2)])]);
% nll = -log(pval);
% nll = sum(nll) / ndata;
    
    
    