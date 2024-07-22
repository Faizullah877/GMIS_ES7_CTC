function report = organize_and_get_report_JPEG_AI(input)

report = struct();
set_info = input.set_info;
no_of_images = set_info.no_of_images;
ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;
q_v_str = sprintf('%03d', q_v);
ai_compressed_folder = fullfile(input.data_folder, "bit");
ai_decompressed_folder = fullfile(input.data_folder, "rec");

for imgNo = 1: no_of_images
    image_name_wo_ext = set_info.images_details.NameWoExt(imgNo);
    no_rows = set_info.images_details.ImageHeight(imgNo);
    no_cols = set_info.images_details.ImageWidth(imgNo);
    compressed_fname = sprintf("VM_%s_%s.bits", image_name_wo_ext, q_v_str);
    source_compressed_path = fullfile(ai_compressed_folder, compressed_fname);
   
    cmp_file_exists = isfile(source_compressed_path);
    if cmp_file_exists
        copyfile(source_compressed_path, input.compressed_folder);
        jpg_image_file_info=dir(source_compressed_path);
        compressed_img_size = jpg_image_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));
    else
        compressed_img_size = NaN;
        img_bpp = NaN;
    end
    img_enc_time  = NaN; 
    decompressed_fname = sprintf("VM_%s_%s.png", image_name_wo_ext, q_v_str);
    source_decomp_path = fullfile(ai_decompressed_folder, decompressed_fname);
    dec_file_exists = isfile(source_decomp_path);
    if dec_file_exists
        decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name_wo_ext));
        ffmpeg_exe = fullfile(input.codec_folder, "ffmpeg.exe");
        conv_cmd=sprintf('%s -i %s %s',ffmpeg_exe, source_decomp_path, decoded_image_ppm);
        [status, ~] = system(conv_cmd);
        if(status)
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", source_decomp_path, decoded_image_ppm);
            error('conversion failed');
        end
    end
    img_dec_time = NaN; %toc(tDstartimg);
    % get bpp and encoding time
    ImagesResult(imgNo,:) ={imgNo, image_name_wo_ext, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end

% replace missing bpp with average
mean_bpp = mean(ImagesResult.bpp, 'omitnan');
nanIndices = isnan(ImagesResult.bpp);
ImagesResult.bpp(nanIndices) = mean_bpp;

% replace misssing filesize with average
mean_size = mean(ImagesResult.CompressedImageSize, 'omitnan');
nanIndices = isnan(ImagesResult.CompressedImageSize);
ImagesResult.CompressedImageSize(nanIndices) = mean_size;

report.set_enc_time = sum(ImagesResult.("enc_time"));
compressed_set_size = sum(ImagesResult.("CompressedImageSize"));
report.compressed_set_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
report.set_dec_time = sum(ImagesResult.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult;

end


