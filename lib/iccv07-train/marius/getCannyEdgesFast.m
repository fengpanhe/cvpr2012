%
%   [edgeComponents, edgeImage, I_orig, BW_orig] = getCannyEdges(I, params)
%
%    Returns the most significant edges. These are edges that are more likely to
%    belong to object boundaries.
%
%   Input:
%     I: Input image array
%     params: parameter structure from read_params
%
%   Output:
%     edgeComponents:
%     %edgeImage:
%     I_orig:
%     BW_orig:
%
%   Exceptions:
%     None
%

function [edgeComponents, I_orig, BW_orig] = getCannyEdgesFast(I, params)

t1 = cputime;
disp('Starting get edges -----');

%% Parameters ------------------------------------------------------------

%USE_BERKLEY = 0;

USE_BERKLEY = get_num_param(params, 'USE_BERKLEY')

%USE_COLOR = 0;

USE_COLOR = get_num_param(params, 'USE_COLOR')

%maxSizeEdges = 320;

maxSizeEdges = get_num_param(params, 'maxSizeEdges');

%KD_radius = 5;

KD_radius = get_num_param(params, 'KD_radius');

%angle_sigma = 10;

angle_sigma = get_num_param(params, 'angle_sigma');

angle_thresh = 3*angle_sigma;

%minRatio1 = 0.01;

minRatio1 = get_num_param(params, 'minRatio1');

%lambdaRatio = 0.5;

lambdaRatio = get_num_param(params, 'lambdaRatio');

%minBW = 0.01;

minBW = get_num_param(params, 'minBW');

%useKDtree = 0;

useKDtree = get_num_param(params, 'useKDtree');

%------------------------------------------------------------------------

%% Preprocessing

I_orig = I;

I  = double(I);

if size(I,3) == 3

    I_aux = (I(:,:,1) + I(:,:,2) + I(:,:,3))/3;

else

    I_aux = I;

end

[nRows, nCols, aux] = size(I);

m = max([nRows, nCols]);

if m > maxSizeEdges

    I = imresize(I, maxSizeEdges/m);
    I_orig = imresize(I_orig, maxSizeEdges/m);
    I_aux = imresize(I_aux, maxSizeEdges/m);

end

[nRows, nCols, aux] = size(I);

t = cputime;
BW = edge(I_aux,'canny', [0.01, 0.3]);
fprintf('CANNY %f', cputime-t);

%figure(21),imshow(BW);
BW_orig = BW;

%pack;

BW = log(BW./minBW);

%compute gradient ....
t = cputime;
[Ix, Iy] = gradient(conv2(I_aux, fspecial('gaussian', 7, 1.5), 'same'));
grad_angle = (atan2(Iy, Ix)) * 180/pi;
fprintf('GRADIENT %f', cputime-t);

ind_edge = find(BW > 0);

x_edge = floor((ind_edge - 1)/nRows) + 1;

y_edge = mod(ind_edge - 1, nRows) + 1;

nEdges = length(ind_edge);

mapEdges = zeros(nRows, nCols);

mapEdges(ind_edge) = 1:nEdges;

max_entries = 100000;

spc_i_indices   = zeros(max_entries, 1);
spc_j_indices   = zeros(max_entries, 1);
spc_affinities  = zeros(max_entries, 1);

components =  1:nEdges;

n_nz = 0;

%use KD trees for edges ------------------------

disp(' Storing edges into a KD tree ');

if useKDtree
    [tmp, tmp, TreeRoot] = kdtree( [y_edge, x_edge], []);
else
    idx = zeros(nRows,nCols);
    idx((x_edge-1)*nRows + y_edge) = [1:nEdges];
end

%% connect edges based on smoothness and build matrix M

disp(' building the matrix ');

t = cputime;
nEdges

for i = 1:nEdges

    if useKDtree
        [locs, dist, index] = kdrangequery( TreeRoot, [y_edge(i), x_edge(i)], KD_radius);
    else
%        [locs, dist, index] = getNeighbors2([y_edge, x_edge],[y_edge(i), x_edge(i)],KD_radius);
        [locs, dist, index] = getNeighbors3([y_edge, x_edge],[y_edge(i), x_edge(i)],idx,KD_radius);
    end

    nNeighbors = length(index);

    for j = 1:nNeighbors

        if index(j) <= i

            continue;

        end

        diff1 = abs(grad_angle(y_edge(i), x_edge(i)) - grad_angle(locs(j,1), locs(j,2)));

        angle_diff1 = min(360 - diff1, diff1);

        angle_i_j = atan2(locs(j,1) - y_edge(i), locs(j,2) - x_edge(i)) * 180/pi;

        diff2 = abs(grad_angle(y_edge(i), x_edge(i)) - angle_i_j);

        angle_diff2 = abs(90 - min(360 - diff2, diff2));

        angle_j_i = atan2(y_edge(i) - locs(j,1), x_edge(i) - locs(j,2)) * 180/pi;

        diff3 = abs(grad_angle(locs(j,1), locs(j,2)) - angle_j_i);

        angle_diff3 = abs(90 - min(360 - diff3, diff3));

        angle_diff4 = abs(min(360 - diff2, diff2) -  min(360 - diff3, diff3));

        if angle_diff1 < angle_thresh && angle_diff2 < angle_thresh &&  angle_diff3 < angle_thresh % ...
            % && angle_diff4 < angle_thresh

            n_nz = n_nz + 1;

            spc_i_indices(n_nz)  = i;
            spc_j_indices(n_nz)  = index(j);

            a1 = 4.5 - ((angle_diff1^2)/(2*angle_sigma^2));

            a2 = 4.5 - ((angle_diff2^2)/(2*angle_sigma^2));

            a3 = 4.5 - ((angle_diff3^2)/(2*angle_sigma^2));

            a4 = 0; %4.5 - ((angle_diff4^2)/(2*angle_sigma^2));

            spc_affinities(n_nz) = 5*(BW(y_edge(i), x_edge(i))+BW(locs(j,1), locs(j,2)))+(a1 + a2 + a3 + a4);

            %connect these edges into components

            c_i = getComponent(components, mapEdges(y_edge(i), x_edge(i)));

            c_j = getComponent(components, mapEdges(locs(j,1), locs(j,2)));

            if c_i > c_j

                components(c_i) = c_j;

            else

                components(c_j) = c_i;

            end


        end

    end

end

fprintf('FIRST LOOP %f', cputime-t);

%finish up the connected components stuff

t = cputime;
for i=1:nEdges
    components(i) = getComponent(components, i);
end

fprintf('COMPONENTS %f', cputime-t);

t = cputime;

n_nz = n_nz + 1;

spc_i_indices(n_nz)  = nEdges;
spc_j_indices(n_nz)  = nEdges;
spc_affinities(n_nz) = 0;

n_nz

spc_i_indices = spc_i_indices(1:n_nz);
spc_j_indices = spc_j_indices(1:n_nz);
spc_affinities = spc_affinities(1:n_nz);

t = cputime;
M = spconvert([spc_i_indices spc_j_indices spc_affinities]);

M = M + M';
fprintf('MATRIX  %f', cputime-t);
%% find significant curves

n = size(M,1);

v = ones(n,1);

maxIter = 50;

for i = 1:maxIter

    v = (M*v);

end

lambda = (M*v) ./ v;

lambda(find(isnan(lambda))) = 0;

u_comp = unique(components);

nNodes = ones(n,1);

for  i = 1:length(u_comp)

    f = find(components == u_comp(i));

    nNodes(f) = length(f);

end

m = max(lambda);

lambda(find(lambda < m*lambdaRatio)) = 0;

scores = lambda .* nNodes;

m = max(scores);

f = find(scores >= m*minRatio1);

edgeImage = zeros(nRows, nCols,3);
edgeComponents = zeros(nRows,nCols);

edgeComponents(ind_edge(f)) = components(f);

 u_comp = unique(components(f));
 
 nColors = length(u_comp);
 
 c = rand(nColors, 3);
 
 lambda_color = zeros(max(u_comp),1);
 lambda_color(u_comp) = 1:nColors;
 
 colors = c(lambda_color(components(f)),:);
 
 temp = zeros(nRows, nCols);
 
 temp(ind_edge(f)) = colors(:,1);
 edgeImage(:,:,1) = temp;
 
 temp(ind_edge(f)) = colors(:,2);
 edgeImage(:,:,2) = temp;
 
 temp(ind_edge(f)) = colors(:,3);
 edgeImage(:,:,3) = temp;

fprintf('REST %f', cputime-t);

fprintf('TOTAL %f', cputime-t1);

return

%% auxiliary function used for the connected components algorithm

function i = getComponent(comp, i)

while comp(i) ~= i

    i = comp(i);

end

return