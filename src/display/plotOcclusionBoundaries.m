function plotOcclusionBoundaries(bndinfo, blabels)
    % im = drawOcclusionBoundaries(im, bndinfo, blabels)
    % blabels(i) = (0, 1, or 2) for off, left, or right
    %

    if numel(blabels) == bndinfo.ne * 2
        blabels = blabels(1:end / 2) + 2 * blabels(end / 2 + 1:end);
    end

    hold on
    imsize = bndinfo.imsize;

    color = 'k';

    indices = bndinfo.edges.indices;

    arrowdist = ceil(sqrt(imsize(1).^2 + imsize(2).^2) / 10);

    for k = 1:numel(indices)

        if blabels(k) > 0

            ind = double(indices{k});
            [ey, ex] = ind2sub(bndinfo.imsize(1:2), ind);
            npix = numel(ind);

            narrows = ceil(npix / arrowdist);

            epos = ceil((1:narrows) / (narrows + 1) * npix);

            for j = 1:numel(epos)

                [ay, ax] = ind2sub(imsize(1:2), ind(epos(j)));

                if blabels(k) == 1
                    [y1, x1] = ind2sub(imsize(1:2), ind(max(epos(j) - 10, 1)));
                    [y2, x2] = ind2sub(imsize(1:2), ind(min(epos(j), npix)));
                else % blabels(k)==2;
                    [y1, x1] = ind2sub(imsize(1:2), ind(min(epos(j) + 10, npix)));
                    [y2, x2] = ind2sub(imsize(1:2), ind(min(epos(j), npix)));
                end

                theta = atan2(y2 - y1, x2 - x1);
                %             if blabels(k)==2
                %                 theta = mod(theta+pi, 2*pi);
                %             end

                asx = ax - cos(theta);
                asy = ay - sin(theta);

                %[ax,ay] = dsxy2figxy(gca, ax, ay);
                %[asx, asy] = dsxy2figxy(gca, asx, asy);
                ax = (ax - 1) / (imsize(2) - 1); 
                ay = 1 - (ay - 1) / (imsize(1) - 1);
                asx = (asx - 1) / (imsize(2) - 1); 
                asy = 1 - (asy - 1) / (imsize(1) - 1);
                ax = max(ax, 0);
                asx = max(asx, 0);

                % plot(ex, ey, 'Color', 1 - color, 'LineWidth', 3);
                plot(ex, ey, 'Color', color, 'LineWidth', 1);
                annotation('arrow', [asx ax], [asy ay], 'Color', color,'LineStyle', 'none', ...
                'HeadStyle', 'vback3' ,'HeadWidth', 17, 'HeadLength', 10);
            end

        end

    end

end
