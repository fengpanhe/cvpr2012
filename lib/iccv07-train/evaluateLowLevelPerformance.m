function result = evaluateLowLevelPerformance(y, pB, w)

ne = numel(y)/2;

if ~exist('w', 'var') || isempty(w) 
    w = ones(ne, 1);
end
w = w / sum(w);
    

valid = (y(1:ne)>=0) & (y(ne+1:end)>=0);

y = (y(1:ne)>0) + 2*(y(ne+1:end)>0);


%%  Edge on/off pr

ty = (y(valid) > 0);

[sval, sind] = sort(1-pB(valid, 1), 'descend');
sy = ty(sind);
fpcsum = cumsum((sy==0).*w(valid));
tpcsum = cumsum((sy==1).*w(valid));
pr = tpcsum ./ (tpcsum + fpcsum);
result.pr.p = pr;
result.pr.r = tpcsum / tpcsum(end);
result.conf = sval;
result.onOffLabel = sy;
result.roc.tp = tpcsum / tpcsum(end);
result.roc.fp = fpcsum / fpcsum(end);
result.numOn = tpcsum(end);
result.numOff = fpcsum(end);


%% Figure/ground accuracy

valid = valid & (y>0);

ty = y(valid);

ty_not = 1*(ty==2) + 2*(ty==1);

rightind = find(valid) + ty*ne;
wrongind = find(valid) + ty_not*ne;

result.labels_fg = ty;

result.fgaccuracy = sum((pB(rightind) > pB(wrongind)) .* w(valid));

% count a tie as half correct
result.fgaccuracy = result.fgaccuracy + ...
    sum((pB(rightind)==pB(wrongind)).*w(valid))*0.5;

result.fgaccuracy = result.fgaccuracy / sum(w(valid));