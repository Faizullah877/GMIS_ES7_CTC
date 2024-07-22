clc
clear all

date_time_str = datestr(datetime('now', 'Format','yyyy-MM-dd HH:mm:ss.SSS'));
date_time_str = erase(date_time_str, ["-", ":", " "]);

diary_name = strcat("diary_", date_time_str,".txt");
diary(diary_name);
diary on

%% Set Experiment parameters
exp_params = struct();
exp_params.set_name = "GolfClubHead_CityScape_10thF";
exp_params.images_files_ext = "*.tiff"; % * must be provided.  
exp_params.set_path = fullfile("D:\GMIS_DataSet\MutliFocusV2", exp_params.set_name); 
exp_params.experiment_folder = "D:\GMIS_EXPs\Evaluate_10Fv2";
exp_params.results_mat_fname = 'Metrics.mat';
exp_params.QA_Metrics = ["PSNR_Y", "PSNR_YUV", "PSNR_HVS_M_Y", "MSSSIM_Y"];

% Anchors_parameters Cell formats
% 1st cell = Anchor Name
% 2nd cell = Quality params array
% 3rd cell = Raw File format for Anchor, a directory with the cell value will be created in [experiment_folder/set_folder/ground_truth] folder
% 4th cell = Is Raw in required format exists or not. if not then the software will generate with same name as 3 cell inside the experiment_folder/set_name/GroundTruth/3Parameter. 
exp_params.Anchors_parameters = {... 
    % {'JPEG_AI',                           [012, 025, 050, 075, 100],          "PNG", true}...
    {'SJU_Arch_JPEG_AI',                  [012, 025, 050, 075, 100],          "PNG", true}...
    % {'SJU_Arch_JPEG1_DPCMPIXELS_RDOPT',   [0.15, 0.30, 0.50, 0.8, 1.0],       "PPM", true, 5} ...
    % {'JPEG1',                             [7, 20, 40, 60, 80, 85, 90, 92, 95], "PPM", false}... %[7, 20, 40, 60, 80, 85]
    % {'JPEG1_Arithmetic',                  [7, 20, 40, 60, 80, 85, 90, 92, 95], "PPM", true}... %[7, 20, 40, 60, 80, 85]
    % {'JPEG2000',                          [7, 20, 40, 60, 80, 85, 90, 92, 95], "PPM", true}... %[7, 20, 40, 60, 80, 85]
    % {'JPEGXL',                            [7, 20, 40, 60, 80, 85, 90, 92, 95], "PPM", true}... %[10, 30, 50, 70, 85, 90]
    % {'SJU_Arch_DPCM_QDCT',                [5, 10, 25, 50, 75, 80, 85, 90, 95], "JPEG-YCbCr", false}... %[10, 30, 50, 70, 85]
    % {'SJU_Arch_DPCM_PIXELSv2',            [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}...}
    % {'SJU_Arch_JPEG2000v2',               [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}... %[10, 25, 50, 75, 90, 95, 5, 93]
    % {'SJU_Arch_JPEG_XLv2',                [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}... %[10, 25, 50, 75, 90, 95, 5, 93]
    % {'VVC_VTM_Intra',                     [50,  40, 32, 25, 20, 16],          "yuv420p10le_images", false}... % [50, 35, 25, 20, 15, 13, 11, 9]
    % {'VVC_VVenC_Inter',                   [50, 35, 30, 25, 20, 15, 10],       "yuv420p10le_video", true}... % [50, 35, 30, 25, 20, 15, 10, 5]
    % {'SJU_Arch_DPCM_PIXELSv2',            [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}...}
    % {'SJU_Arch_JPEG2000v2',               [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}... %[10, 25, 50, 75, 90, 95, 5, 93]
    % {'SJU_Arch_JPEG_XLv2',                [5, 10, 20, 40, 60, 75, 85, 90],    "PPM", true}... %[10, 25, 50, 75, 90, 95, 5, 93]
    % {'VVC_VVenC_Inter_1_thread',          [50, 35, 30, 25],                   "yuv420p10le",    false}... %[50, 35, 30, 25, 20, 15, 10, 5],
    % {'VVC_VVenC_Inter_8_thread',          [50, 35, 30, 25],                   "yuv420p10le", true}... %[50, 35, 30, 25, 20, 15, 10, 5],
    % {'VVC_VVenC_Intra',                   [50, 35, 30, 25, 20, 15, 10],       "yuv420p10le_images", false}... % [50, 35, 30, 25, 20, 15, 10, 5]
    % {'VVC_VTM_Inter',                       [50, 35, 25, 20, 15, 12, 10, 8],    "yuv420p10le_video", false}... % [50, 35, 25, 20, 15, 12, 10, 8]
    % {'H264_AVC', [25, 30, 35], "PPM", true}...
    % {'H264_AVCintra', [25, 30, 35], "PPM", true}...
    % {'SJU_Arch_DPCM_PIXELS', [10, 25, 50, 75, 90, 95, 5, 18, 40, 70, 80, 85], "JPEG-YCbCr", true}...}
    };

exp_params.do_discard = false;

%% Add paths to functions etc
addpath(genpath('./functions'));
addpath(genpath('./GMIS_IQA'));
addpath(genpath('./UNSW_v4_functions'));

exp_params.dir_struct = prepare_directories(exp_params);


if exist(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'file')
    load(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname));
else
    M=struct();
end

if ~exist(fullfile(exp_params.dir_struct.set_data_folder, 'dataset_info.mat'), 'file')
    dataset_info = get_dataset_info(exp_params.set_path, exp_params.set_name, exp_params.images_files_ext);
    save(fullfile(exp_params.dir_struct.set_data_folder, 'dataset_info.mat'), 'dataset_info');
else
    load(fullfile(exp_params.dir_struct.set_data_folder, 'dataset_info.mat'));
end
exp_params.set_info = dataset_info;
M.set_info = dataset_info;
%% loop through each anchor
no_of_anchors = numel(exp_params.Anchors_parameters);

for ancNo = 1:no_of_anchors
    anch_P = exp_params.Anchors_parameters{ancNo}; 
    
    if ~anch_P{4}
        data_format = anch_P{3};
        exp_params = convert_dataset_fileformat(exp_params, data_format);
    end

    M = compression_pipeline_v2(exp_params, anch_P, M);
    save(fullfile(exp_params.dir_struct.results_folders, exp_params.results_mat_fname), 'M');
end

if exp_params.do_discard
    rmdir(exp_params.dir_struct.set_data_folder, "s");
end
diary off

%% create the necessary directories. 
function dir_struct = prepare_directories(exp_params)
    
    current_dir = pwd;
    dir_struct.set_path = exp_params.set_path;
    dir_struct.codec_folder = fullfile(current_dir, "codecs");
    dir_struct.gmis_iqa = fullfile(current_dir, "GMIS_IQA");
    
    dir_struct.experiment_folder = exp_params.experiment_folder;
    if(~exist(dir_struct.experiment_folder, 'dir')), mkdir(dir_struct.experiment_folder); end

    dir_struct.results_folders = fullfile(dir_struct.experiment_folder, "Results");
    if(~exist(dir_struct.results_folders, 'dir')), mkdir(dir_struct.results_folders); end
      
    dir_struct.set_data_folder = fullfile(dir_struct.experiment_folder, exp_params.set_name);
    if(~exist(dir_struct.set_data_folder, 'dir')), mkdir(dir_struct.set_data_folder); end
    
    dir_struct.ground_truth = fullfile(dir_struct.set_data_folder, "GroundTruth");
    if(~exist(dir_struct.ground_truth, 'dir')), mkdir(dir_struct.ground_truth); end
end