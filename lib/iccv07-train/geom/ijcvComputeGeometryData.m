function [imdata, spdata] = ijcvComputeGeometryData(im, imsegs)

spdata = mcmcGetSuperpixelData(im, imsegs); 
imdata = mcmcComputeImageData(im, imsegs);
