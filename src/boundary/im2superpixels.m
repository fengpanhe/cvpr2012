function imsegs = im2superpixels(im)

    fn1 = [tempname '.ppm'];
    fn2 = [tempname '.ppm'];
    segcmd = './bin/segment 0.8 100 100';
    
    imwrite(im, fn1);
    system([segcmd ' ' fn1 ' ' fn2]);
    imsegs = processSuperpixelImage(fn2);
    
    
    delete(fn1);
    delete(fn2);
end