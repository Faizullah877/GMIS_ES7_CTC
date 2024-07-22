function report = encode_decode_VVC_VTM(input)

report = struct();
set_info = input.set_info;
set_name = input.set_name;
q_v = input.q_value;
no_of_images = set_info.no_of_images;

encoder_exe = fullfile(input.codec_folder, "VVC_VTM", "EncoderApp.exe");
ffmpeg_exe = fullfile(input.codec_folder, "ffmpeg.exe");

ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};

if strcmp(input.config, "intra")
    cfg_file = fullfile(input.codec_folder, "VVC_VTM","cfg", "encoder_intra_vtm.cfg");
    for imgNo = 1 : no_of_images
        no_rows = set_info.images_details.ImageHeight(imgNo);
        no_cols = set_info.images_details.ImageWidth(imgNo);
        image_name = set_info.images_details.NameWoExt(imgNo);
        yuv_file_name = sprintf("%s_%dx%d_1Hz_prec10_P420.yuv", image_name, no_cols, no_rows);
        input_yuv_image = fullfile(input.input_raw, yuv_file_name); % file naming format is followed in conversion.
        if ~exist(input_yuv_image, 'file')
            error("Input File does not exists");
        end
        bin_image_file = fullfile(input.compressed_folder, sprintf("%s.bin", image_name));
        decoded_image_yuv = fullfile(input.decompressed_folder, sprintf("Rec_%s.yuv", image_name));
        enc_cmd =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f 1 -fr 1 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 -o %s', encoder_exe, cfg_file, input_yuv_image, no_cols, no_rows, bin_image_file,  q_v, decoded_image_yuv);
        tCstartimg= tic;
        [status, outE] = system(enc_cmd);
        tCcomplete = toc(tCstartimg);
        img_enc_time = (tCcomplete/4)*3;
        img_dec_time = tCcomplete/4;
        if(status)
            fprintf("VVC VTM => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_image);
            disp(outE);
            error('compression failed');
        end
        compressed_file_info=dir(bin_image_file);
        compressed_img_size = compressed_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));
        ppm_output=fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name));
        conv_cmd = sprintf('%s -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', ffmpeg_exe, no_cols, no_rows, decoded_image_yuv, ppm_output);
        [status, ~] = system(conv_cmd);
        if(status)
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_image_yuv, ppm_output);
            error('conversion failed');
        end
        delete(decoded_image_yuv);
        ImagesResult(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
    end
    report.set_enc_time = sum(ImagesResult.("enc_time"));
    compressed_set_size = sum(ImagesResult.("CompressedImageSize"));
    report.compressed_set_size = compressed_set_size;
    report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
    report.set_dec_time = sum(ImagesResult.("dec_time"));
    report.set_dec_folder = input.decompressed_folder;
    report.images_result_table = ImagesResult;
    
else
    if ~set_info.same_resolution
        fprintf("The codec {VVC-VTM-Inter} does not support sets {%s} with images having different resolution.", set_name);
        error("Compression Failed");
    end
    no_rows = set_info.images_details.ImageHeight(1);
    no_cols = set_info.images_details.ImageWidth(1);
    input_yuv_video = fullfile(input.input_raw, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", set_name, no_cols,no_rows, no_of_images)); % file naming format is followed in conversion.
    bin_video_file = fullfile(input.compressed_folder, sprintf("%s.bin", set_name));
    decoded_video_yuv = fullfile(input.decompressed_folder, sprintf("Rec_%s.yuv", set_name));
    cfg_file = fullfile(input.codec_folder, "VVC_VTM","cfg", "encoder_randomaccess_vtm.cfg");
    enc_cmd =sprintf('%s -c %s -i %s -wdt %d -hgt %d -b %s -f %d -fr 30 -q %d --InputBitDepth=10 --InputChromaFormat=420 --ChromaFormatIDC=420 --TemporalSubsampleRatio=1 --ConformanceWindowMode=1 --MinSearchWindow=1 -o %s', encoder_exe, cfg_file, input_yuv_video, no_cols, no_rows, bin_video_file, no_of_images, q_v, decoded_video_yuv);
    tCstartimg= tic;
    [status, outE] = system(enc_cmd);
    report.set_enc_time = toc(tCstartimg);
    if(status)
        fprintf("VVC VTM => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
        error('compression failed');
    end

    compressed_file_info=dir(bin_video_file);
    report.compressed_file = bin_video_file;
    report.compressed_set_size = compressed_file_info.bytes;
    report.set_bpp = (compressed_file_info.bytes*8) / (double(no_cols)*double(no_rows)*double(no_of_images));
    num_of_bits_Table = get_images_bpp_vvc_out(outE);
    time_taken = extract_process_time(outE);
    report.set_enc_time = time_taken/4*3;
    report.set_dec_time = time_taken/4;

    ppm_output=fullfile(input.decompressed_folder,"dec_image_%03d.ppm");
    conv_cmd = sprintf('%s -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', ffmpeg_exe, no_cols, no_rows, decoded_video_yuv, ppm_output);
    [status, ~] = system(conv_cmd);
    if(status)
        fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_video_yuv, decoded_rgb_files);
        error('conversion failed');
    end
    delete(decoded_video_yuv);
    for imgNo = 1: no_of_images
        image_name = set_info.images_details.NameWoExt(imgNo);
        row = num_of_bits_Table.Image_No == imgNo;
        num_of_bits = num_of_bits_Table.Num_of_Bits(row);
        num_of_pixels = no_rows*no_cols;
        compressed_img_size = num_of_bits/8; % to convert to bytes from bits
        img_bpp = num_of_bits/num_of_pixels;
        img_enc_time = report.set_enc_time/no_of_images;
        img_dec_time = report.set_dec_time/no_of_images;
        movefile(fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm",imgNo)), fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name)));
        ImagesResult(imgNo,:) ={imgNo, image_name,  uint16(no_cols), uint16(no_rows),  compressed_img_size, img_bpp, img_enc_time, img_dec_time};
    end
    report.set_dec_folder = input.decompressed_folder;
    report.images_result_table = ImagesResult;
end
end

function time = extract_process_time(charVector)
% Split the char vector into lines
lines = strsplit(charVector, '\n');


for i = 1:length(lines)
    line = strtrim(lines{i});  % Remove leading and trailing white spaces
    if startsWith(line, 'Total Time:')  % Check if the line starts with 'Total Time:'
        tokens = strsplit(line);  % Split the line into words
        time = str2double(tokens{3});  % Add the number after 'Total Time:' to time
    end
end
end

