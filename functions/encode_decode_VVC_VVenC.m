function report = encode_decode_VVC_VVenC(input)

        report = struct();
        set_name = input.set_name;
        q_v = input.q_value;
        no_cols = input.width;
        no_rows = input.height;
        no_of_images = input.no_of_images;
        no_of_threads = input.no_of_threads;
        
        encoder_exe = fullfile(input.codec_folder, "VVC_VVenC", "vvencFFapp.exe");

        input_yuv_video = fullfile(input.input_raw, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", set_name, no_cols,no_rows, no_of_images)); % file naming format is followed in conversion.
        % input_yuv_video = input.input_raw;
        bin_video_file = fullfile(input.compressed_folder, sprintf("%s.266", set_name));
        rnd_cfg_file = fullfile(input.codec_folder, "VVC_VVenC/cfg/randomaccess_slow.cfg");
        seq_cfg_file = fullfile(input.codec_folder, "VVC_VVenC/cfg/sequence.cfg");

        if(strcmp(input.config, "inter"))
            enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --MinSearchWindow 1 --framerate 30 -t %d", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_video, no_cols, no_rows, no_of_images, bin_video_file, no_of_threads);
        else  
            enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --framerate 30 --IntraPeriod=1  -t %d", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_video, no_cols, no_rows, no_of_images, bin_video_file, no_of_threads);            
        end


        tCstartimg= tic;
        [status, outE] = system(enc_cmd);
        report.set_enc_time = toc(tCstartimg);
        if(status)
            fprintf("VVC VVC_VVenC => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
            error('compression failed');
        end

        compressed_file_info=dir(bin_video_file);
        report.compressed_file = bin_video_file;
        report.compressed_file_size = compressed_file_info.bytes;
        report.set_bpp = (compressed_file_info.bytes*8) / (input.width*input.height*input.no_of_images);

        num_of_bits_Table = get_images_bpp_vvc_out(outE);
        time_taken = extract_enc_time(outE);
        report.set_enc_time = time_taken;
        % decoding
        decoder_exe = fullfile(input.codec_folder, "VVC_VVenC/vvdecapp.exe");
        decoded_video_yuv = fullfile(input.decompressed_folder, sprintf("Rec_%s.yuv", set_name));

        dec_cmd = sprintf("%s -b %s -o %s ", decoder_exe,  bin_video_file, decoded_video_yuv);

        tDstartimg= tic;
        [status, ~] = system(dec_cmd);
        report.set_dec_time = toc(tDstartimg);
        if(status)
            fprintf("VVC vvdecapp => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
            error('compression failed');
        end

        % tConStartT = tic;
        ffmpeg_exe = fullfile(input.codec_folder, "ffmpeg.exe");
        ppm_output=fullfile(input.decompressed_folder,"dec_image_%03d.ppm");
        conv_cmd = sprintf('%s -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', ffmpeg_exe, no_cols, no_rows, decoded_video_yuv, ppm_output);
        [status, ~] = system(conv_cmd);
        if(status)
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_video_yuv, decoded_rgb_files);
            error('conversion failed');
        end
        % tx = toc(tConStartT);
        % fprintf("\t\tConversion Done in %f sec.\n", tx);
        delete(decoded_video_yuv);

        report.set_dec_folder = input.decompressed_folder;


        ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
        ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
        for imgNo = 1: no_of_images
            row = num_of_bits_Table.Image_No == imgNo;
            num_of_bits = num_of_bits_Table.Num_of_Bits(row);
            num_of_pixels = no_rows*no_cols;
            compressed_img_size = num_of_bits/8; % to convert to bytes from bits
            img_bpp = num_of_bits/num_of_pixels;
            img_enc_time = report.set_enc_time/no_of_images;
            img_dec_time = report.set_dec_time/no_of_images;
            image_name = sprintf("image%03d", imgNo);
            ImagesResult_T1(imgNo,:) ={imgNo, image_name,  uint16(no_cols), uint16(no_rows),  compressed_img_size, img_bpp, img_enc_time, img_dec_time};
        end
        
        report.images_result_table = ImagesResult_T1;
end
function time = extract_enc_time(charVector)
% Split the char vector into lines
lines = strsplit(charVector, '\n');


for i = 1:length(lines)
    line = strtrim(lines{i});  % Remove leading and trailing white spaces
    if startsWith(line, 'vvencFFapp [info]: Total Time:')  % Check if the line starts with 'Total Time:'
        tokens = strsplit(line);  % Split the line into words
        time = str2double(tokens{5});  % Add the number after 'Total Time:' to time
    end
end
end

