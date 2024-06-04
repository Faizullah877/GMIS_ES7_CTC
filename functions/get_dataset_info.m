function set_details = get_dataset_info(set_path, set_name,  filesExts)


images_list = dir(fullfile(set_path, filesExts));
images_list=images_list(~ismember({images_list.name},{'.','..'}));
n = length(images_list);

img_names = strings(n,1);
img_size = zeros(n,1);
img_width = zeros(n,1);
img_height = zeros(n,1);
img_ncolors = zeros(n,1);
img_ftype = strings(n,1);
for i = 1: n
    if ~isfolder(images_list(i).name)
        img_names(i,1) = images_list(i).name;
        img_size(i,1) = images_list(i).bytes;
        img = imread(fullfile(images_list(i).folder, images_list(i).name));
        img_height(i,1) = size(img,1);
        img_width(i,1) = size(img,2);
        img_ncolors(i,1) = size(img,3);
        [~, ~, ext] = fileparts(fullfile(images_list(i).folder, images_list(i).name));
        img_ftype(i,1) = ext;
    end
end
img_details = table(img_names, img_size, img_width, img_height, img_ncolors, img_ftype);
set_details.images_details = img_details;
set_details.no_of_images = n;
set_details.total_size = sum(img_details.img_size);


if all(img_details.img_width==img_details.img_width(1))
    set_details.same_width = 1;
else
    set_details.same_width = 0;
end
if all(img_details.img_height==img_details.img_height(1))
    set_details.same_height = 1;
else
    set_details.same_height = 0;
end

if  (set_details.same_width == 1)&&(set_details.same_height == 1)
    set_details.same_resolution = 1;
else
    set_details.same_resolution = 0;
end


set_details.set_name = set_name;
set_details.set_path = set_path;
end

