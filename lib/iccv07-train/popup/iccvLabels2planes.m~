function [gplanes, vplanes, gim, vim] = ...
    iccvLabels2planes(bndinfo, labels, hy, im)
% [gplanes, vplanes, gim, vim] = ...
%    APPlabels2planes(vlabels, hlabels, hy, segimage, image)
% Converts a label map into a set of planes
%
% planes(nplanes).{inds(npoints), vertices(1:2, npoints), label}
% gim and vim are the maps of which pixels are vertical and ground

ALL_BLOCKS = 1;

regions = labels.regions;
rgeom = labels.regionGeometry;
echain = labels.edgeletChain;
spLR = labels.chainLR;  
ordering = labels.depthOrder;
alljuncts = [bndinfo.edges.junctions ; bndinfo.edges.junctions(:, [2 1])];

%% Estimate fit to ground for each vertical object
% For each vertical region, determine how it is connected with ground
% Describe ground-boundary with a polyline
[sval, sind] = sort(ordering, 'ascend');
vind = sind(rgeom(sind)==2); % indices of vertical regions from front to back
for k = 1:numel(vind)
    eind = find((spLR(:, 1)==vind(k)) | (spLR(:, 2)==vind(k)));        
     
    % get all junctions for each vertical object and record the other
    % region adjacent to each
    nj = 1000;
    juncts = zeros(nj,1);
    jregion = zeros(nj,1);
    nj = 0;
    for ec = eind'
        for e = echain{ec}
            nj = nj+1;
            juncts(nj) = alljuncts(e,1);
            jregion(nj) = spLR(ec, spLR(ec,:)~=k);
        end
    end
    juncts = juncts(nj);
    jregion = jregion(nj);
    
    % find contact points (initially assuming object is whole)
    
    if any(jisgnd)
        jpos = bndinfo.junctions.position(juncts, :);
        

%% Warp vertical objects to frontal view
% Save separate texture maps for each vertical object.


vim = zeros(size(segimage));
gim = zeros(size(segimage));
gplanes = cell(0, 1);
vplanes = cell(0, 1); 

labelmap = zeros(size(segimage));
[height, width] = size(labelmap);
nsegs = length(vlabels);

for h = 1:height
    for w = 1:width
        s = segimage(h, w);
        if strcmp(vlabels{s}, '000')
            labelmap(h, w) = 1;
        elseif strcmp(vlabels{s}, '090') | strcmp(vlabels{s}, '045')
            labelmap(h, w) = 2;
        end
    end
end


% get all vertical planes
pcount = 0;
labels = {'000', '090'};

for y = 2:length(labels)    

    [L, nplanes] = get_object_map(labelmap==y);

    %[L, nplanes] = bwlabel(labelmap==y, 4);
            
    %nplanes = max(L(:));
    
    for m = 1:nplanes               

        inds = find(L==m);        
        %disp(length(inds))
        if numel(inds) > 2
                        
            yinds = mod(inds-1, height)+1;
            xinds = floor((inds-1)/height) + 1;            
            % get information about whole plane
                   
            gnd_pts = APPcreateGroundPoints(L==m, labelmap==1);                                                       
            %disp(length(gnd_pts))
            
            lfparam = [];
            %disp(size(gnd_pts, 1))
            
            if size(gnd_pts, 1) > 4 % enough for line fit
                %lfparam = piecewise_linear_spline2(gnd_pts(:, 2), gnd_pts(:, 1), 2);    
                lfparam = APPfitGroundHough(gnd_pts(:, 2), gnd_pts(:, 1), [height width]);   
                
            else
                labelmap(inds) = (yinds < (1-hy));
                inds = [];
            end
                

            if size(gnd_pts, 1) > 0

                for s = 1:size(lfparam, 1)
					pcount = pcount + 1;
                
                    %vplanes(pcount).gfit = lfparam(s, :);                            
                    vplanes(pcount).gfit = lfparam;
                    vplanes(pcount).label = labels{y};                    
                    %vplanes(pcount).disconnect = lfparam(s, 6);
                    
                    if s==size(lfparam, 1)
                        x1 = lfparam(end, 3)*(width-1)+1;
                        x2 = width;
                    elseif s==1
                        x1 = 1;
                        x2 = lfparam(1, 4)*(width-1)+1;
                    else
                        x1 = lfparam(s, 3)*(width-1)+1;
                        x2 = lfparam(s, 4)*(width-1)+1;
                    end
                    if s < size(lfparam, 1) && ~lfparam(s, 6)
                        x2 = min(x2 + 1, width);
                    end
                    
                    ti = find(xinds>=x1 & xinds<=x2);                    
                    sxinds = xinds(ti);
                    syinds = yinds(ti);
                    
                    if s==1 
                        ti2 = find(xinds<(lfparam(1, 3)*(width-1)+1));
                        %vim(inds(ti2)) = 1;
                        xinds2 = xinds(ti2);
                        yinds2 = yinds(ti2);
                        m = lfparam(1, 1);
                        b = lfparam(1, 2);
                        ti3 = find(yinds2 <= round((m*(xinds2-1)/(width-1) + b)*(height-1)+1));
                        vim(inds(ti2(ti3))) = 1;                        
                    end
                    if s==size(lfparam, 1)
                        ti2=  find(xinds>lfparam(end, 4)*(width-1)+1);
                        %vim(inds(ti2)) = 1;
                        xinds2 = xinds(ti2);
                        yinds2 = yinds(ti2);
                        m = lfparam(end, 1);
                        b = lfparam(end, 2);
                        ti3 = find(yinds2 <= round((m*(xinds2-1)/(width-1) + b)*(height-1)+1));
                        vim(inds(ti2(ti3))) = 1;                        
                    end
                    if x1~=x2
     	                y1 = min(yinds(ti));                
                        tx1 = lfparam(s, 3)*(width-1)+1;
                        tx2 = lfparam(s, 4)*(width-1)+1;
                        y21 = (lfparam(s, 3)*lfparam(s, 1) + lfparam(s, 2))*(height-1)+1;
 		                y22 = (lfparam(s, 4)*lfparam(s, 1) + lfparam(s, 2))*(height-1)+1;                     
                        %disp([num2str(s) ': ' num2str(round([y21 y22]))])
                        for tx = round(tx1):round(tx2)
                            ty2 = min(max(y21+(tx-x1)/(x2-x1)*(y22-y21), 1), height); 
                            vim(round(y1):round(ty2), tx) = 1;
                        end
                    %else
                        %vim(ti) = 1;
                    end
                    miny = min(syinds);
					maxy = max(syinds);
                    minx = min(sxinds);
					maxx = max(sxinds);  
                    if x1~=x2
                        maxy = max([y21 y22 maxy]);
                    end
                    
                    count = 0;
                    step = 15;
                    vplanes(pcount).label = '090';
                    vplanes(pcount).sinfo = cell(ceil((maxy-miny)/step)*ceil((maxx-minx)/step), 1);
                    for j = maxy-step:-step:miny
                        ti = find((syinds >= j) & (syinds <= j+step));
                        %minx = min(sxinds(ti));
                    	%maxx = max(sxinds(ti));        
                        %disp(maxx-minx)                        
                        for i = minx:step:maxx-step
                            if ALL_BLOCKS || any((sxinds(ti)>=i) & (sxinds(ti)<=i+step)) %new
                                count = count + 1;
                                vplanes(pcount).sinfo{count}.x = [i       i+step  i+step i];
                                vplanes(pcount).sinfo{count}.y = [j+step  j+step  j      j];                           
                            end
                        end
                        i = max([minx i]);
                        if i+step < maxx
                            if ALL_BLOCKS || any((sxinds(ti)>=i) & (sxinds(ti)<=maxx)) %new
                                count = count + 1;
                                vplanes(pcount).sinfo{count}.x = [i       maxx    maxx   i];
                                if length(vplanes(pcount).sinfo{count}.x) < 4
                                    disp(['error: ' num2str(i) ' - ' num2str(maxx)])
                                end
                                vplanes(pcount).sinfo{count}.y = [j+step  j+step  j      j];                          
                            end
                        end
                    end  
                    %disp(length(ti))
                    if isempty(j) | j ~= miny 
                        j = min([j maxy-step]);
                        ti = find((syinds <= j) & (syinds >= miny));
                        %minx = min(sxinds(ti));
                    	%maxx = max(sxinds(ti));          
                        for i = minx:step:maxx-step
                            if ALL_BLOCKS || any((sxinds(ti)>=i+step) & (sxinds(ti)<=maxx)) %new
                                count = count + 1;
                                vplanes(pcount).sinfo{count}.x = [i    i+step  i+step  i   ];
                                vplanes(pcount).sinfo{count}.y = [j    j       miny    miny];                         
                            end
                        end
                        i = max([minx i]);
                        if ALL_BLOCKS || any((sxinds(ti)>=i) & (sxinds(ti)<=maxx)) %new
                            count = count + 1;                       
                            vplanes(pcount).sinfo{count}.x = [i  maxx  maxx  i   ];
                            vplanes(pcount).sinfo{count}.y = [j+step  j+step     miny  miny];                            
                        end
                    end
                    vplanes(pcount).sinfo = vplanes(pcount).sinfo(1:count); 
                end
                if count == 0
                    vplanes(pcount) = [];
                    pcount = pcount - 1;
                end                
                
            end            
                     
        end
    end
end


% get single ground plane
inds = find((labelmap~=0) & ~vim);
np = 0;
if length(inds > 0)
	yinds = mod(inds-1, height)+1;
	xinds = floor((inds-1)/height) + 1;

	minx = min(xinds);
	maxx = max(xinds);
	miny = min(yinds);
	maxy = max(yinds);

    gplanes(1).label = '000';    

   
    step = 7;
    count = 0;
    gplanes(1).sinfo = [];
    for j = maxy-step:-step:miny               
        ti = find((yinds <= j) & (yinds >= miny));        
        for i = minx:step:maxx-step                                    
            if ALL_BLOCKS || any((xinds(ti)>=i) & (xinds(ti)<=i+step)) %new
                count = count + 1;
                gplanes(1).sinfo{count}.x = [i       i+step  i+step i];
                gplanes(1).sinfo{count}.y = [j+step  j+step  j      j];                           
            end
        end
        i = max([minx i]);
        if i+step < maxx
            if ALL_BLOCKS || any((xinds(ti)>=i) & (xinds(ti)<=maxx)) %new
                count = count + 1;
                gplanes(1).sinfo{count}.x = [i       maxx    maxx   i];
                if length(gplanes(1).sinfo{count}.x) < 4
                    disp(['error: ' num2str(i) ' - ' num2str(maxx)])
                end
                gplanes(1).sinfo{count}.y = [j+step  j+step  j      j];                          
            end
        end
    end  
    %disp(length(ti))
    if isempty(j) | j ~= miny 
        j = min([j maxy-step]);
        ti = find((yinds <= j) & (yinds >= miny));
        %minx = min(sxinds(ti));
    	%maxx = max(sxinds(ti));          
        for i = minx:step:maxx-step
            if ALL_BLOCKS || any((xinds(ti)>=i) & (xinds(ti)<=i+step)) %new
                count = count + 1;
                gplanes(1).sinfo{count}.x = [i    i+step  i+step  i   ];
                gplanes(1).sinfo{count}.y = [j    j       miny    miny];                         
            end
        end
        i = max([minx i]);
        if ALL_BLOCKS || any((xinds(ti)>=i) & (xinds(ti)<=maxx)) %new
            count = count + 1;                       
            gplanes(1).sinfo{count}.x = [i  maxx  maxx  i   ];
            gplanes(1).sinfo{count}.y = [j+step  j+step     miny  miny];                            
        end
    end            
    
    gplanes(1).sinfo = gplanes(1).sinfo(1:count);

else
    gplanes = [];
end

vim = vim & (labelmap~=0);
gim = (labelmap~=0) & (~vim);
%figure(4), imshow(gim)
gim = imdilate(gim, strel('square',20));


% end labels2planes_s10



function [L, nplanes] = get_object_map(ymap)
L1 = imerode(ymap, strel('square', 11));
L1 = imerode(ymap, strel('square', 11));
L1 = imerode(ymap, strel('square', 11));
[L1, nplanes1] = bwlabel(L1, 4);
%figure(5), imshow(label2rgb(L1, 'jet'))
L1 = imdilate(L1, strel('square', 15));
L1 = imdilate(L1, strel('square', 15));    
L1 = imdilate(L1, strel('square', 15));     
L1 = L1.*(ymap);


L = L1;
nplanes = nplanes1;


   
    
    
