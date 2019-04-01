function test_im2boundariesTopLevel(dataset_path, image_name)

    image_file = fullfile(dataset_path, 'images', strcat(image_name, '.jpg'));

    im = imread(image_file);
    [bndinfo, pbim, gconf, bndinfo_all] = im2boundariesTopLevel(double(im) / 255);
    bndinfo_file = fullfile('result', strcat(image_name, '_ex.mat'));
    save(bndinfo_file, 'bndinfo', 'pbim', 'gconf', 'bndinfo_all');

    dest_dir_path = fullfile('result');
    writeContourDepthResults(bndinfo_file, image_file, dest_dir_path, 1);
end
