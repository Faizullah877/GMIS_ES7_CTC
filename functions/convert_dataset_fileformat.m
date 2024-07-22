function exp_params = convert_dataset_fileformat(exp_params, data_format)

tCstart = tic;
set_info = exp_params.set_info;
ffmpeg_exe = fullfile(exp_params.dir_struct.codec_folder, "ffmpeg.exe");
output_folder = fullfile(exp_params.dir_struct.ground_truth, data_format);
exp_params.ground_truth_folder = output_folder;
no_of_images =set_info.no_of_images;
if(~exist(output_folder, 'dir')), mkdir(output_folder); end
if strcmp(data_format, "yuv444p10le_video") || strcmp(data_format, "yuv420p10le_video")
%% for video files
    % to make sure the order of images in video file, files are copied and
    % renamed to a sequencial order.
    tmp_folder = fullfile(exp_params.dir_struct.ground_truth, "temp");
    if(~exist(tmp_folder, 'dir')), mkdir(tmp_folder); end
    no_of_images = exp_params.set_info.no_of_images;
    for i = 1:no_of_images
        srcFile = fullfile(set_info.set_path, set_info.images_details.NameWdExt(i));
        ext = set_info.images_details.Ext(i);
        destFile = fullfile(tmp_folder, sprintf("image_%03d%s",  i, ext));
        copyfile(srcFile, destFile);
    end
    input = fullfile(tmp_folder, sprintf("image_%%03d%s", ext));
    we = set_info.images_details.ImageWidth(1);
    he = set_info.images_details.ImageHeight(1);
    use_ffmpeg_cmd = true;
    if strcmp(data_format , "yuv444p10le_video")
        output_file = fullfile(output_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P444.yuv", exp_params.set_name, we, he, no_of_images));
        cmd = sprintf("%s -hide_banner -i %s -pix_fmt yuv444p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", ffmpeg_exe, input, output_file);
    elseif strcmp(data_format , "yuv420p10le_video")
        output_file = fullfile(output_folder, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", exp_params.set_name, we,he, no_of_images));
        cmd = sprintf("%s -hide_banner -i %s -pix_fmt yuv420p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", ffmpeg_exe, input, output_file);
    else
        error("Unable to convert to video file for Pixel Format %s", data_format);
    end

    if use_ffmpeg_cmd
        [status, out] = system(cmd);
        if status
            disp(out);
            error("failed to convert the set");
        end
        rmdir(tmp_folder, "s");
    end
else
    %% for images files
    for imgNo = 1 : no_of_images
        image_name = set_info.images_details.NameWoExt(imgNo);
        image_name_wd_ext = set_info.images_details.NameWdExt(imgNo);
        input_file = fullfile(set_info.set_path, image_name_wd_ext);
        we = set_info.images_details.ImageWidth(imgNo);
        he = set_info.images_details.ImageHeight(imgNo);
        use_ffmpeg_cmd = true;
        if strcmp(data_format, "PPM")
            output_file = fullfile(output_folder, sprintf("%s.ppm", image_name));
            cmd = sprintf("%s -hide_banner -i %s -y %s", ffmpeg_exe, input_file, output_file);
        elseif strcmp(data_format , "yuv444p10le_images")
            output_file = fullfile(output_folder, sprintf("%s_%dx%d_1Hz_prec10_P444.yuv", image_name, we, he));
            cmd = sprintf("%s -hide_banner -i %s -pix_fmt yuv444p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", ffmpeg_exe, input_file, output_file);
        elseif strcmp(data_format , "yuv420p10le_images")
            output_file = fullfile(output_folder, sprintf("%s_%dx%d_1Hz_prec10_P420.yuv", image_name, we, he));
            cmd = sprintf("%s -hide_banner -i %s -pix_fmt yuv420p10le -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709  -color_trc bt709 -colorspace bt709 -y %s", ffmpeg_exe, input_file, output_file);
        elseif strcmp(data_format, "JPEG-YCbCr")
            rgb_img = imread(fullfile(exp_params.set_path, image_name_wd_ext));
            [he, we, ~] = size(rgb_img);
            ycbcr_img = rgb2jpegycbcr(rgb_img);
            sju_imwrite_yuv(ycbcr_img, fullfile(output_folder, sprintf("%s_%dx%d_1Hz_prec8_P444.yuv", image_name, we,he)));
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
    end

end

conversion_time = toc(tCstart);
fprintf(['   Conversion done in  ' num2str(conversion_time) ' sec']);
exp_params.ground_truth_folder = output_folder;
end


