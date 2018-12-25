function createBoundaryInfo(fn, imdir, pbdir, savedir)

for f = 1:numel(fn)
    
    disp(num2str(f))        

    savename = [savedir '/' strtok(fn{f}, '.') '_seg.mat'];
    
    if exist(savename, 'file') 
        disp(['Warning: overwriting file ' savename ' because it already exists.']);
    end
        
    if 1
        im = im2double(imread([imdir fn{f}]));    
        load([pbdir strtok(fn{f}, '.') '_pb.mat']);        

        wseg = pb2wseg(pball, 100000);
        
        [edges, juncts, neighbors, wseg, avecol, polyfragments, polyparams] = ...
            seg2fragments(double(wseg), im, 25);
        bndinfo = processBoundaryInfo3(wseg, edges, neighbors);

        bndinfo.imname = fn{f};

        for k = 1:numel(polyparams)
            bndinfo.edges.polyparams{k} = single(polyparams{k});
        end
    
        disp(['Saving ' savename ' ...'])
        
        imind = f;
        save(savename, 'bndinfo', 'imind');
    end
    
end


    