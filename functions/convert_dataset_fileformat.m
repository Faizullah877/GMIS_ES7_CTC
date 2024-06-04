function exp_params = convert_dataset_fileformat(exp_params, data_format)

tCstart = tic;
dir_struct = exp_params.dir_struct;
image_name_template = char(exp_params.set_info.images_details.img_names(1));

dotPosition = strfind(image_name_template, '.');
numDigits = 0; 
for i = dotPosition-1:-1:1
    if isstrprop(image_name_template(i), 'digit') % check if the character is a digit
        numDigits = numDigits + 1; % increment the count
    else
        break;
    end
end

image_name_base = string(image_name_template(1:i));
ext = string(image_name_template(dotPosition:end));
input = sprintf("%s/%s%%0%dd%s" , exp_params.set_path,image_name_base, numDigits,ext);

output_folder = fullfile(exp_params.dir_struct.ground_truth, data_format);
if(~exist(output_folder, 'dir')), mkdir(output_folder); end
we = exp_params.set_info.images_details.img_width(1);
he = exp_params.set_info.images_details.img_height(1);
no_of_images = exp_params.set_info.no_of_images;
use_ffmpeg_cmd = true;
if strcmp(data_format, "PPM")
    output_file = sprintf("%s/%s%%0%dd.ppm", output_folder, image_name_base, numDigits);
    cmd = sprintf("ffmpeg -hide_banner -i %s -y %s", input, output_file);
elseif strcmp(data_format , "yuv444p10le")
    output_file = fullfile(output_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P444.yuv", exp_params.set_name, we, he, no_of_images));
    cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt yuv444p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input, output_file);
elseif strcmp(data_format , "yuv420p10le")
    output_file = fullfile(output_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", exp_params.set_name, we,he, no_of_images));
    cmd = sprintf("ffmpeg -hide_banner -i %s -pix_fmt yuv420p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", input, output_file);
elseif strcmp(data_format, "JPEG-YCbCr")
    for img=1:no_of_images
        rgb_img = imread(fullfile(exp_params.set_path, exp_params.set_info.images_details.img_names(img)));
        [he, we, ~] = size(rgb_img);
        ycbcr_img = rgb2jpegycbcr(rgb_img);
        sju_imwrite_yuv(ycbcr_img, fullfile(output_folder, sprintf("image_%03d_%dx%d_1Hz_prec8_P444.yuv", img, we,he)));
    end
    use_ffmpeg_cmd = false;
else
    error("Unable to convert to video file for Pixel Format %s", data_format);
end

if use_ffmpeg_cmd
    [status, ~] = system(cmd);
    if status
        error("failed to convert the set");
    end
end
conversion_time = toc(tCstart);
fprintf(['   Conversion done in  ' num2str(conversion_time) ' sec']);

end


