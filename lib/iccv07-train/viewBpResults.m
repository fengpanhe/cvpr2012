function viewBpResults(bndinfo, Y, pB, pC, adjlist)

ne = bndinfo.ne;

[factors, f2var] = getContourPotentials(pB, pC, adjlist, bndinfo, 1E-3);
nnodes = 3*ones(ne, 1);
[mllab, bel] = maxBeliefPropBethe(factors, f2var, nnodes, 0.025, 0.15, Inf);
bel = cell2num(bel')';  

pcim = zeros(size(bndinfo.wseg));   
for k = 1:size(pB,1)/2
    pcim(bndinfo.edges.indices{k}) = 1-pB(k, 1);
end
figure(1), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray    


pcim = zeros(size(bndinfo.wseg));   
for k = 1:size(pB,1)/2
    pcim(bndinfo.edges.indices{k}) = 1-bel(k,1);% 1-pB(k, 1);
end
figure(2), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray

pcim = zeros(size(bndinfo.wseg));   
for k = 1:size(pB,1)/2
    pcim(bndinfo.edges.indices{k}) = (1-bel(k,1))-(1-pB(k,1));% 1-pB(k, 1);
end
figure(5), hold off, imagesc(ordfilt2(pcim,9,ones(3))), axis image, colormap gray


%% Edge on/off pr
y = (Y(1:ne)>0) | (Y(ne+1:end)>0);

valid = (Y(1:ne)>=0) & (Y(ne+1:end)>=0);
y = y(valid);

[sval, sind] = sort(1-pB(valid, 1), 'descend');
sy = y(sind);
fpcsum = cumsum(sy==0);
tpcsum = cumsum(sy==1);
pr = tpcsum ./ (tpcsum + fpcsum);
%figure(3), hold off, plot(fpcsum / fpcsum(end), tpcsum / tpcsum(end), 'b');
figure(3), hold off, plot(tpcsum / tpcsum(end), pr, 'b');


[sval, sind] = sort(1-bel(valid, 1), 'descend');
sy = y(sind);
fpcsum = cumsum(sy==0);
tpcsum = cumsum(sy==1);
pr = tpcsum ./ (tpcsum + fpcsum);
%figure(3), hold on, plot(fpcsum / fpcsum(end), tpcsum / tpcsum(end), 'g'); hold off
figure(3), hold on, plot(tpcsum / tpcsum(end), pr, 'g');
axis([0 1 0 1])
%keyboard

figure(4), plot(sval, tpcsum / tpcsum(end))

%% Edge direction pr
y = (Y(1:ne)>0) + 2*(Y(ne+1:end)>0);
valid = (Y(1:ne)>=0) & (Y(ne+1:end)>=0) & (y>0);

ny = 1*(y==2) + 2*(y==1);

tmppB = [pB(1:ne, :) pB(ne+1:end, 2)];
accuracy(1) = mean(tmppB(find(valid) + y(valid)*ne) > tmppB(find(valid) + ny(valid)*ne));
accuracy(2) = mean(bel(find(valid) + y(valid)*ne) > bel(find(valid) + ny(valid)*ne));

disp(num2str(accuracy));



