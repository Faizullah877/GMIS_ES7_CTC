function set_details = get_dataset_info(set_path, set_name,  filesExts)

set_details = struct();

images_list = dir(fullfile(set_path, filesExts));
if numel(images_list)==0
    error("Unable to get the files in input directory for conversion");
end
images_list=images_list(~ismember({images_list.name},{'.','..'}));
no_of_images = length(images_list);
ImagesDetails = table('Size', [no_of_images 9], 'VariableTypes', {'uint16', 'string', 'string','string', 'uint64', 'uint16', 'uint16', 'uint64','uint16'});
ImagesDetails.Properties.VariableNames = {'ImageNo', 'NameWoExt','NameWdExt', 'Ext', 'FileSize', 'ImageWidth', 'ImageHeight', 'NumPixels','NumColors'};

for imgNo = 1: no_of_images
    if ~isfolder(images_list(imgNo).name)
        image_file_name = images_list(imgNo).name;
        [~, name_wo_ext, ext] = fileparts(image_file_name);
        file_size = images_list(imgNo).bytes;
        img = imread(fullfile(images_list(imgNo).folder, images_list(imgNo).name));        
        [no_rows, no_cols, no_colors] = size(img);
        num_pixels = no_cols * no_rows;
        ImagesDetails(imgNo,:) ={imgNo, name_wo_ext, image_file_name, ext, file_size, no_cols, no_rows, num_pixels, no_colors};
    end
end
set_details.images_details = ImagesDetails;
set_details.no_of_images = no_of_images;
set_details.total_size = sum(ImagesDetails.FileSize);
set_details.num_pixels = sum(ImagesDetails.NumPixels);

if all(ImagesDetails.ImageWidth==ImagesDetails.ImageWidth(1))
    set_details.same_width = 1;
else
    set_details.same_width = 0;
end
if all(ImagesDetails.ImageHeight==ImagesDetails.ImageHeight(1))
    set_details.same_height = 1;
else
    set_details.same_height = 0;
end

if  (set_details.same_width == 1)&&(set_details.same_height == 1)
    set_details.same_resolution = 1;
else
    set_details.same_resolution = 0;
end

set_details.file_ext = ImagesDetails.Ext(1);
set_details.set_name = set_name;
set_details.set_path = set_path;
end

