function imsegs = im2superpixels(im)

prefix = num2str(floor(rand(1)*10000000));
prefix = '1';
fn1 = ['./result/tmp/tmpim' prefix '.ppm'];
fn2 = ['./result/tmp/tmpimsp' prefix '.ppm'];
segcmd = './bin/segment 0.8 100 100';

imwrite(im, fn1);
system([segcmd ' ' fn1 ' ' fn2]);
imsegs = processSuperpixelImage(fn2);


delete(fn1);
delete(fn2);