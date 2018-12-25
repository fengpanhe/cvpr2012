function [pr, roc, fgacc] = summarizeLowLevelPerformance(result)

%% average PR curve

allr = (0:0.005:1);
p = zeros(size(allr));
for f = 1:numel(result)        
    
    c = 0;
    for r = allr
        
        c = c + 1;
        ind = find(result(f).pr.r==r);
        if isempty(ind)                
            [tmp, ind] = min(abs(result(f).pr.r - r));
        else
            ind = ind(1);
        end
        p(c) = p(c) + result(f).pr.p(ind)/numel(result);
        
    end
end
pr.p = p;
pr.r = allr;


%% average ROC curve

allfp = (0:0.005:1);
tp = zeros(size(allfp));
for f = 1:numel(result)        
    
    c = 0;
    for fp = allfp
        
        c = c + 1;
        ind = find(result(f).roc.fp==fp);
        if isempty(ind)                
            [tmp, ind] = min(abs(result(f).roc.fp - fp));
        else
            ind = ind(1);
        end
        tp(c) = tp(c) + result(f).roc.tp(ind)/numel(result);
        
    end
end
roc.fp = allfp;
roc.tp = tp;


%% Average figure ground accuracy

fgacc = mean([result(:).fgaccuracy]);


    