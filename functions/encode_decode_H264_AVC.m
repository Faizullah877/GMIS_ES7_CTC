function report = encode_decode_H264_AVC(input)

        report = struct();

        encoder_exe = fullfile(input.codec_folder, "ffmpeg.exe");

        orig_files = dir(fullfile(input.input_raw, "*.ppm"));
        [~, image_name, ~] = fileparts(orig_files(1).name);
        
        ext_str = sprintf("%%03d%s", ".ppm");

        enc_input =fullfile(input.input_raw, erase(image_name,"001")+ ext_str );
        compressed_file = fullfile(input.compressed_folder, sprintf('%s_q%d.264',input.set_name,input.q_value ));

        if (strcmp(input.config , "inter"))
            enc_command = sprintf('%s -y -s %dx%d -r %d -i %s -bf 1 -c:v libx264 -crf %d  %s',encoder_exe, input.width, input.height, 1, enc_input, input.q_value, compressed_file); % -y is to overwrite if output file already exist.
        else
            enc_command = sprintf('%s -y -s %dx%d -r %d -i %s -g 1 -c:v libx264 -crf %d  %s',encoder_exe, input.width, input.height, 1, enc_input, input.q_value, compressed_file); % -y is to overwrite if output file already exist.
        end
        tCstartimg= tic;

        [status, ~] = system(enc_command);
        report.set_enc_time = toc(tCstartimg);
        if status
            error("Error: Compression Failed.");
        end

        mp4_file_info=dir(compressed_file);
        report.compressed_file = compressed_file;
        report.compressed_file_size = mp4_file_info.bytes;
        report.set_bpp = (mp4_file_info.bytes*8) / (input.width*input.height*input.no_of_images);

        % decoding
        decoder_exe = encoder_exe;
        ppm_output=fullfile(input.decompressed_folder,"dec_image_%03d.ppm");
        dec_command = sprintf ('%s -y -i %s %s', decoder_exe, compressed_file, ppm_output );
        tDstartimg = tic;
        [status, ~] = system(dec_command);
        report.set_dec_time = toc(tDstartimg);
        if status
            error("Error: De-compression Failed.");
        end
        report.set_dec_folder = input.decompressed_folder;

end

