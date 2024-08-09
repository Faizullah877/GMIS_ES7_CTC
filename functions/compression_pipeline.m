function M = compression_pipeline(exp_params, codec_params, M)
    %COMPRESSION_PIPELINE Summary of this function goes here
    %   Detailed explanation goes here

    %general parameters
    set_info = exp_params.set_info;
    set_name = exp_params.set_name;
    QA_metrics = exp_params.QA_Metrics;
    do_discard = exp_params.do_discard;   
    codec_exp_time = tic;    
    codec = codec_params{1};
    q_list = codec_params{2};
    data_format = codec_params{3};   
    no_of_QAmetrics = numel(QA_metrics);    
    quality_count = numel(q_list);
    dir_struct = exp_params.dir_struct;  
    no_of_images = set_info.no_of_images;    
    set_codec_results = struct();
    ImagesNames = cellstr(set_info.images_details.NameWoExt');

    % codec data folder
    codec_data_folder = fullfile(dir_struct.set_data_folder, codec);
    if(~exist(codec_data_folder, 'dir')), mkdir(codec_data_folder); end    
    


    % check for existing average results table for a set.
    [M, append_row] = get_result_table(M, set_name, ImagesNames, codec, quality_count, QA_metrics);
    set_avg_results = M.(set_name).(codec).SET_AVG_RESULTS;

    %% Iterate over rate/q values
    for q_itr = 1: quality_count
        q_v = q_list(q_itr);
        q_fieldName = ['q_' num2str(q_v)];
        q_fieldName = strrep(q_fieldName, ".", "p"); % replace dot with char 'p' in case of floating point q values. 

        % check and skip if the data for the q value is already there.
        if isfield(M, set_name)
            if isfield(M.(set_name), codec)
                if isfield(M.(set_name).(codec), "q_fields")
                    if isfield(M.(set_name).(codec).q_fields, q_fieldName)
                        str = sprintf("\n\tSet Name : [%s] => Codec : [%s] => Q field : [%s] already exists and the experiment is skipped.\n", set_name, codec, q_fieldName);
                        fprintf(str);
                        continue
                    end
                end
            end
        end
    
        % create the folders for the compressed and decompressed data.
        set_q_folder = fullfile(codec_data_folder, q_fieldName);
        if(~exist(set_q_folder, 'dir')), mkdir(set_q_folder); end
        set_compressed_folder = fullfile(set_q_folder, "Compressed");
        if(~exist(set_compressed_folder, 'dir')), mkdir(set_compressed_folder); end
        set_dec_folder = fullfile(set_q_folder, "Decompressed");
        if(~exist(set_dec_folder, 'dir')), mkdir(set_dec_folder); end

        % necessary variables for codec and results
        ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
        ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};    
        Images_IQA_T = double(zeros([no_of_images no_of_QAmetrics])); 
    
        % input structure contain data for codec
        input = struct();
        input.set_name = set_name;
        input.codec_folder = dir_struct.codec_folder;
        if strcmp(codec, "JPEG_AI") | strcmp(codec, "SJU_Arch_JPEG_AI")
            input.input_raw = exp_params.set_path;
        else
            input.input_raw = fullfile(dir_struct.ground_truth, data_format);
        end
        input.set_q_folder = set_q_folder;
        input.compressed_folder = set_compressed_folder;
        input.decompressed_folder = set_dec_folder;
        input.q_value = q_v;
        input.set_info = set_info;
    
        fprintf("\tStart Time: %s.\n", get_current_time_str());
        fprintf("\tCompression Pipeline running for SET : [%s {No. Images = %d}]   CODEC : [%s] and RATE : [%3.2f] => %d/%d\n",set_name,no_of_images, codec, q_v, q_itr, quality_count);
    
    
        q_start_time = tic();
    
    
    
        if strcmp(codec, "H264_AVC")
            input.config = "inter";
            codec_report = encode_decode_H264_AVC(input);
        elseif strcmp(codec, "H264_AVCintra")
            input.config = "intra";
            codec_report = encode_decode_H264_AVC(input);
        elseif strcmp(codec, "JPEG1")
            input.config = "huff";
            codec_report = encode_decode_JPEG1(input);
        elseif strcmp(codec, "JPEG1_Arithmetic")
            input.config = "arith";
            codec_report = encode_decode_JPEG1(input);
        elseif strcmp(codec, "JPEG2000")
            codec_report = encode_decode_JPEG2000(input);
        elseif strcmp(codec, "JPEGXL")
            codec_report = encode_decode_JPEG_XL(input);
        elseif strcmp(codec, "VVC_VVenC_Inter")
            input.config = "inter";
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Inter_1_thread")
            input.config = "inter";
            input.no_of_threads = 1; 
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Inter_8_thread")
            input.config = "inter";
            input.no_of_threads = 8; 
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VVenC_Intra")
            input.config = "intra";
            codec_report = encode_decode_VVC_VVenC(input);
        elseif strcmp(codec, "VVC_VTM_Inter")
            input.config = "inter";
            codec_report = encode_decode_VVC_VTM(input);
        elseif strcmp(codec, "VVC_VTM_Intra")
            input.config = "intra";
            codec_report = encode_decode_VVC_VTM(input);
        elseif strcmp(codec, "SJU_Arch_DPCM_QDCT")
            input.encoding_scheme = "jpeg1_dpcm_qdct_seq";
            codec_report = encode_decode_DPCM_QDCT(input);
            % codec_report = encode_decode_DPCM_QDCT(input);
        elseif strcmp(codec, "SJU_Arch_DPCM_PIXELSv2")
            input.config = "arith";
            input.codec = "JPEG1";
            codec_report = encode_decode_SJU_ARCHs(input);
            % input.encoding_scheme = "jpeg1_dpcm_pixels_seq_idx";
            % codec_report = encode_decode_SJU_ARCH(input);
        elseif strcmp(codec, "SJU_Arch_JPEG1_DPCMPIXELS_RDOPT")
            input.config = "arith";
            input.codec = "JPEG1";
            input.intra_refresh_rate = codec_params{5};
            codec_report = encode_decode_SJU_ARCH_JPEG1_RDOPT(input);
        elseif strcmp(codec, "JPEG_AI")
            input.data_folder = codec_data_folder;
            codec_report = organize_and_get_report_JPEG_AI(input);
        elseif strcmp(codec, "SJU_Arch_JPEG_AI")
            input.data_folder = codec_data_folder;
            codec_report = organize_and_get_report_SJU_Arch_JPEG_AI(input);
        elseif strcmp(codec, "SJU_Arch_JPEG2000v2")
            input.codec = "JPEG2000";
            codec_report = encode_decode_SJU_ARCHs(input);
        elseif strcmp(codec, "SJU_Arch_JPEG_XLv2")
            input.codec = "JPEG_XL";
            codec_report = encode_decode_SJU_ARCHs(input);
        else
            fprintf("Codec [%s] functionality not included yet", codec);
            error("Failed");
        end
    
        codec_time = toc(q_start_time);
    
        fprintf('\t\t=> Set Compressed and Decompressed in %f sec.\n', codec_time);
        set_bpp = (codec_report.compressed_set_size*8) / set_info.num_pixels;
        eval_set_start_time = tic;
        textprogressbar('        => Evaluating : ');
        for imgNo=1:no_of_images
            image_name_wo_ext = set_info.images_details.NameWoExt(imgNo);
            image_name_wd_ext = set_info.images_details.NameWdExt(imgNo);
            input_raw = fullfile(exp_params.set_info.set_path, image_name_wd_ext);
            Orig_IMG = imread(input_raw);
            [no_rows, no_cols, ~] = size(Orig_IMG);
            decoded_ppm_file = fullfile(set_dec_folder, sprintf('%s.ppm', image_name_wo_ext));
            if strcmp(codec, "SJU_Arch_DPCM_QDCT")
                decoded_ppm_file = fullfile(set_dec_folder, sprintf('dec_image_%03d.ppm', imgNo));
            end
            dec_file_exists = isfile(decoded_ppm_file);
            if dec_file_exists
                Recon_IMG = imread(decoded_ppm_file);
                image_IQA_values = evaluate_images_IQA(Orig_IMG, Recon_IMG, QA_metrics);
                if ~isfield(codec_report, "images_result_table")
                    compressed_img_size = codec_report.compressed_set_size/no_of_images;
                    img_bpp = set_bpp;
                    img_enc_time = codec_report.set_enc_time/no_of_images;
                    img_dec_time = codec_report.set_dec_time/no_of_images;
                    ImagesResult_T1(imgNo,:) ={imgNo, image_name_wo_ext, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
                end
                tmp = struct2cell(image_IQA_values);
                v = [tmp{:}];
                Images_IQA_T(imgNo,:) = v;
            else             
                Images_IQA_T(imgNo,:) = NaN;
            end
            textprogressbar(imgNo/no_of_images*100);
        end
        eval_set_time = toc(eval_set_start_time);
        textprogressbar(['   done in  ' num2str(eval_set_time) ' sec']);
    
        if isfield(codec_report, "images_result_table")
            ImagesResult_T1 = codec_report.images_result_table;
        end
    
        fprintf("\tEnd Time: %s.\n\n", get_current_time_str());
    
        m_table = array2table(Images_IQA_T, 'VariableNames',QA_metrics);

        columnNames = m_table.Properties.VariableNames;
        for i = 1:length(columnNames)
            columnName = columnNames{i};
            if isnumeric(m_table.(columnName))
                meanValue = mean(m_table.(columnName), 'omitnan');
                nanIndices = isnan(m_table.(columnName));
                m_table.(columnName)(nanIndices) = meanValue;
            end
        end

        images_results_table = [ImagesResult_T1, m_table];
        if append_row
            existing_rows_count = height(set_avg_results);
            qp_no = existing_rows_count+1;
        else
            qp_no = q_itr;
        end
        set_avg_results(qp_no, :) = {q_v, codec_report.compressed_set_size, set_bpp, codec_report.set_enc_time, codec_report.set_dec_time, ...
            mean(images_results_table.(QA_metrics(1))), ...
            mean(images_results_table.(QA_metrics(2))), ...
            mean(images_results_table.(QA_metrics(3))), ...
            mean(images_results_table.(QA_metrics(4)))};
    
        set_codec_results.(q_fieldName) = images_results_table;
        M.(set_name).(codec).q_fields.(q_fieldName) = images_results_table;
        M.(set_name).(codec).SET_AVG_RESULTS = sortrows(set_avg_results, "SET_BPP");

        M.(set_name).(codec) = fill_results_tables(M.(set_name).(codec), qp_no, q_v, images_results_table, QA_metrics);


        save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');
    
        % delete directorires and files
        if(do_discard)
            rmdir(set_compressed_folder, "s");
            rmdir(set_dec_folder, "s");
            rmdir(set_q_folder, "s");
        end
    
    end

    M.(set_name).(codec).RATE_VALUES = q_list;
    M.(set_name).(codec).IMAGES_NAMES = ImagesNames;
    codec_time_taken = toc(codec_exp_time);
    fprintf("<strong>  => Experiment Done : => SET = [%s],  => CODEC = [%s]  => TIME TAKEN = [%f sec]</strong>\n\n", set_name, codec, codec_time_taken);
    if(do_discard)
        rmdir(codec_data_folder, "s");
    end
end

function timestr = get_current_time_str()
nowTime = datetime('now', 'Format','dd-MMM-uuuu HH:mm:ss');
timestr = string(nowTime);
end

function  [M, append_row]= get_result_table(M, set_name, ImagesNames, codec, quality_count, QA_metrics)
    nq = numel(QA_metrics);
    if isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "SET_AVG_RESULTS")
        append_row = true;
    else
        nnq = 5 +nq;
        SET_AVG_RESULTS = table('Size', [quality_count nnq], 'VariableTypes', repmat("double", 1, nnq));
        variablenames = [["Qvalue", "SET_CompressedSize","SET_BPP", "SET_ENC_TIME", "SET_DEC_TIME"], QA_metrics];
        SET_AVG_RESULTS.Properties.VariableNames = cellstr(variablenames);
        SET_AVG_RESULTS{:,:} = NaN;
        M.(set_name).(codec).SET_AVG_RESULTS = SET_AVG_RESULTS;
        append_row = false;
    end

    table_header = {'Qvalue', 'MEAN', 'MAX', 'MIN', 'STD', 'VAR'};
    table_header = [table_header, ImagesNames];
    nn = numel(table_header);
    % compressed images sizes table
    if ~(isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "ImagesCompressedSizesT"))
        variable_types = repmat("double", 1, nn);
        ImagesCompressedSizesT = table('Size',[quality_count nn], 'VariableTypes', variable_types); %
        ImagesCompressedSizesT.Properties.VariableNames = table_header;
        ImagesCompressedSizesT{:,:} = NaN;
        M.(set_name).(codec).ImagesCompressedSizesT =ImagesCompressedSizesT;
    end
    
    % images bpps table
    if ~(isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "ImagesBPPsT"))
        variable_types = repmat("double", 1, nn);
        ImagesBPPsT = table('Size',[quality_count nn], 'VariableTypes', variable_types); %
        ImagesBPPsT.Properties.VariableNames = table_header;
        ImagesBPPsT{:,:} = NaN;
        M.(set_name).(codec).ImagesBPPsT =ImagesBPPsT;
    end
    % images EncodingTimes table
    if ~(isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "ImagesEncTimesT"))
        variable_types = repmat("double", 1, nn);
        ImagesEncTimesT = table('Size',[quality_count nn], 'VariableTypes', variable_types); %
        ImagesEncTimesT.Properties.VariableNames = table_header;
        ImagesEncTimesT{:,:} = NaN;
        M.(set_name).(codec).ImagesEncTimesT =ImagesEncTimesT;
    end

    % images EncodingTimes table
    if ~(isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), "ImagesDecTimesT"))
        variable_types = repmat("double", 1, nn);
        ImagesDecTimesT = table('Size',[quality_count nn], 'VariableTypes', variable_types); %
        ImagesDecTimesT.Properties.VariableNames = table_header;
        ImagesDecTimesT{:,:} = NaN;
        M.(set_name).(codec).ImagesDecTimesT =ImagesDecTimesT;
    end
    for k = 1: nq
        metric_name = QA_metrics(k);
        if ~(isfield(M, set_name) && isfield(M.(set_name), codec) &&  isfield(M.(set_name).(codec), metric_name))
            variable_types = repmat("double", 1, nn);
            metric_nameT = table('Size',[quality_count nn], 'VariableTypes', variable_types); %
            metric_nameT.Properties.VariableNames = table_header;
            metric_nameT{:,:} = NaN;
            M.(set_name).(codec).(metric_name) =metric_nameT;
        end
    end
end

function CodecResultsT = fill_results_tables(M_codec, qp_no, q_v, images_results_table, QA_metrics)
        M_codec.ImagesCompressedSizesT(all(ismissing(M_codec.ImagesCompressedSizesT), 2), :) = [];
        rowT = get_images_tables_row(images_results_table, q_v, "CompressedImageSize");  
        M_codec.ImagesCompressedSizesT(qp_no, :) = rowT;
        M_codec.ImagesCompressedSizesT = sortrows(M_codec.ImagesCompressedSizesT, "Qvalue");
        M_codec.ImagesBPPsT(all(ismissing(M_codec.ImagesBPPsT), 2), :) = [];
        rowT = get_images_tables_row(images_results_table, q_v, "bpp");  
        M_codec.ImagesBPPsT(qp_no, :) = rowT;
        M_codec.ImagesBPPsT = sortrows(M_codec.ImagesBPPsT, "Qvalue");
        M_codec.ImagesEncTimesT(all(ismissing(M_codec.ImagesEncTimesT), 2), :) = [];
        rowT = get_images_tables_row(images_results_table, q_v, "enc_time");  
        M_codec.ImagesEncTimesT(qp_no, :) = rowT;
        M_codec.ImagesEncTimesT = sortrows(M_codec.ImagesEncTimesT, "Qvalue");
        M_codec.ImagesDecTimesT(all(ismissing(M_codec.ImagesDecTimesT), 2), :) = [];
        rowT = get_images_tables_row(images_results_table, q_v, "dec_time");  
        M_codec.ImagesDecTimesT(qp_no, :) = rowT;
        M_codec.ImagesDecTimesT = sortrows(M_codec.ImagesDecTimesT, "Qvalue");
        nq = numel(QA_metrics);
        for q = 1: nq
            metric = QA_metrics(q);
            if ismember(metric, images_results_table.Properties.VariableNames)
                M_codec.(metric)(all(ismissing(M_codec.(metric)), 2), :) = [];
                rowT = get_images_tables_row(images_results_table, q_v, metric);
                M_codec.(metric)(qp_no, :) = rowT;
                M_codec.(metric) = sortrows(M_codec.(metric), "Qvalue");
            end
        end
        CodecResultsT= M_codec;
end

function rowT = get_images_tables_row(images_results_table, q_v, value)
    rowT = array2table([double(q_v),...
        mean(images_results_table.(value)),...
        max(images_results_table.(value)), ...
        min(images_results_table.(value)), ...
        std(images_results_table.(value)), ...
        var(images_results_table.(value)), ...
        images_results_table.(value)']);         
end

