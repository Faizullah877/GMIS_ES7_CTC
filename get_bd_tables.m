clc
clear all

addpath(genpath('./bjontegaard'));
M_File_Path = "D:\GMIS_EXPs\MultiExposure_v1\Results";

M_file = fullfile(M_File_Path, "Metrics.mat");
plot_all_codecs = false;

tables_folder = fullfile(M_File_Path, "BD_Tables");
if(~exist(tables_folder, 'dir')), mkdir(tables_folder); end

load(M_file);
set_name = "Set_507";

set_results_present = isfield(M, set_name);
if(set_results_present ~= 1)
    fprintf("\t\tResults for Set = %s does not exist in the Metrics file.\n", set_name);
    return;
end

MaxBPP = 8.0; 

% set_table_folder = fullfile(tables_folder, set_name);
% if(~exist(set_table_folder, 'dir')), mkdir(set_table_folder); end

metrices = [ "PSNR_Y", "PSNR_YUV", "PSNR_HVS_M_Y", "MSSSIM_Y"];
% plot_metrics = ["PSNR_RGB", "PSNR_Y",  "PSNR_HVS_M_Y", "MSSSIM_RGB"];

if(plot_all_codecs)
    plot_codecs = string(fieldnames(M.(set_name)));
else
     plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_QDCT", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2", "JPEG_AI", "SJU_Arch_JPEG_AI"];
    
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_QDCT", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2"];
    
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "SJU_Arch_DPCM_PIXELSv2", "JPEG2000", "SJU_Arch_JPEG2000v2", "JPEGXL", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM", "SJU_Arch_JPEG_AI"];  
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000"];
    % codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VTM_Intra", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG2000v2", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM", "SJU_Arch_JPEG_AI"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VTM_Intra", "SJU_Arch_DPCM_QDCT", "SJU_Arch_JPEG2000v1", "SJU_Arch_JPEG_XLv1", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG2000v2", "SJU_Arch_JPEG_XLv2", "JPEG_AI_VM"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VTM_Intra", "SJU_Arch_DPCM_QDCT", "SJU_Arch_JPEG2000", "SJU_Arch_JPEG_XL", "JPEG_AI_VM"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VVenC_Intra"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic","JPEG2000", "JPEGXL", "VVC_VTM_Intra", "VVC_VVenC_Inter", "JPEG_AI_VM", "SJU_Arch_DPCM_QDCT","SJU_Arch_DPCM_PIXELS"];
    % plot_codecs = ["JPEG2000", "JPEGXL", "VVC_VVenC_Inter", "VVC_VVenC_Intra", "VVC_VTM_Inter", "VVC_VTM_Intra"];
    % plot_codecs = ["JPEG2000", "JPEGXL", "SJU_Arch_DPCM_PIXELSv1", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG2000v1", "SJU_Arch_JPEG2000v2", "SJU_Arch_JPEG_XLv1", "SJU_Arch_JPEG_XLv2"];
    % plot_codecs = ["JPEG1", "JPEG1_Arithmetic", "JPEG2000", "JPEGXL", "SJU_Arch_DPCM_PIXELSv2", "SJU_Arch_JPEG1_DPCMPIXELS_RDOPT"];
end
codec_names = { ...
    % Codec Name ,            , codec name
    {"JPEG1"                 ,  "JPEG 1"               } ...
    {"JPEG1_Arithmetic"      ,  "JPEG 1 Arithmethic"   } ...
    {"JPEG2000"              ,  "JPEG 2000"            } ...
    {"JPEGXL"                ,  "JPEG XL"              } ...
    {"VVC_VVenC_Inter"       ,  "VVC inter"            } ...
    {"VVC_VVenC_Intra"       ,  "VVC intra"            } ...
    {"VVC_VTM_Inter"         ,  "VVC inter"            } ...
    {"VVC_VTM_Intra"         ,  "VVC intra"            } ...
    {"SJU_Arch_DPCM_PIXELSv1",  "SJU-Arch JPEG 1"      } ...
    {"SJU_Arch_DPCM_PIXELSv2",  "SJU-Arch JPEG 1 DPCM-PIXELS"      } ...
    {"SJU_Arch_DPCM_QDCT"    ,  "SJU-Arch JPEG 1 DPCM-QDCT"      } ...
    {"SJU_Arch_JPEG2000v1"   ,  "SJU-Arch JPEG 2000"   } ...
    {"SJU_Arch_JPEG2000v2"   ,  "SJU-Arch JPEG 2000"   } ...
    {"SJU_Arch_JPEG_XLv1"    ,  "SJU-Arch JPEG XL"     } ...
    {"SJU_Arch_JPEG_XLv2"    ,  "SJU-Arch JPEG XL"     } ...
    {"JPEG_AI"            ,  "JPEG AI"              } ...
    {"SJU_Arch_JPEG_AI"      ,  "SJU-Arch JPEG AI"     } ...
    {"SJU_Arch_JPEG1_DPCMPIXELS_RDOPT",  "SJU-Arch JPEG1 Rdopt"}
    };
xls_filename = fullfile(tables_folder, sprintf("%s.xlsx", set_name));
bd_rate_sheet = "BD_rates";
bd_quality_sheet = "BD_quality";

for cd =1:length(metrices)
    metric = metrices(cd);
    [BD_rateT, BD_snrT] = get_bjontegaard_tables(M.(set_name), plot_codecs, codec_names, metric, MaxBPP);
    [no_rows, no_cols] = size(BD_rateT);
    C1 = 'A';
    R1 = ((cd-1)*(no_rows+3))+2;
    C2 = char('A' + no_cols +1);
    R2 = R1 + no_rows+2;
    Title_Range = sprintf("%c%d", C1, R1-1);
    writematrix(metric, xls_filename, 'Sheet',bd_rate_sheet, 'Range', Title_Range);
    writematrix(metric, xls_filename, 'Sheet',bd_quality_sheet, 'Range', Title_Range);
    Range = sprintf("%c%d:%c%d",C1, R1, C2,R2);
    writetable(BD_rateT, xls_filename, 'Sheet',bd_rate_sheet, 'Range',Range, 'WriteVariableNames',true,'WriteRowNames',true);
    writetable(BD_snrT, xls_filename, 'Sheet',bd_quality_sheet, 'Range',Range, 'WriteVariableNames',true,'WriteRowNames',true);
    % file_name = fullfile(set_table_folder, sprintf("BD_Rate_%s.txt", metric));
    % writetable(BD_rateT, file_name);
    % file_name = fullfile(set_table_folder, sprintf("BD_Quality_%s.txt", metric));
    % writetable(BD_snrT, file_name);
end


function [BD_rateT, BD_snrT] = get_bjontegaard_tables(M_set, codecs_list, codec_names, metric_name, MaxBPP)
no_of_codecs = numel(codecs_list);
var_names = strings(1, no_of_codecs);
for a=1:no_of_codecs
    codec_name = codecs_list(a);
    matches = cellfun(@(x) strcmp(x{1}, codec_name), codec_names);
    codec_names1 = codec_names(matches);
    codec_names1 = codec_names1{1};
    var_names(a) = codec_names1{2};
end
v_type = "double";
variable_types = repmat(v_type, 1, no_of_codecs);
BD_rateT = table('Size', [no_of_codecs, no_of_codecs], 'VariableTypes', variable_types,  'VariableNames', var_names, 'RowNames', var_names);
BD_snrT = table('Size', [no_of_codecs, no_of_codecs], 'VariableTypes', variable_types,  'VariableNames', var_names, 'RowNames', var_names);

for i=1:no_of_codecs
    for j=i+1:no_of_codecs
        codec1 = codecs_list(i);
        bpp1 = M_set.(codec1).SET_AVG_RESULTS.SET_BPP;
        bpp1 = bpp1(bpp1<MaxBPP);
        no_ele = length(bpp1);
        dist1 = M_set.(codec1).SET_AVG_RESULTS.(metric_name);
        dist1 = dist1(1:no_ele);
        if strcmp(metric_name, "MSSSIM_Y")
            dist1 = -10 *log10(1-dist1);
        end

        codec2 = codecs_list(j);
        bpp2 = M_set.(codec2).SET_AVG_RESULTS.SET_BPP;
        bpp2 = bpp2(bpp2<MaxBPP);
        no_ele = length(bpp2);
        dist2 = M_set.(codec2).SET_AVG_RESULTS.(metric_name);
        dist2 = dist2(1:no_ele);
        if strcmp(metric_name, "MSSSIM_Y")
            dist2 = -10 *log10(1-dist2);
        end
        % [bd_rateA, IoU] = bd_pchip(bpp1, dist1, bpp2, dist2, false);
        drate = bjontegaard(bpp1, dist1, bpp2, dist2, 'rate');
        dsnr = bjontegaard(bpp1, dist1, bpp2, dist2, 'dsnr');
        BD_rateT{var_names(i), var_names(j)} = drate;

        BD_snrT{var_names(i), var_names(j)} = dsnr;
    end
end
% BD_rateT(ismember(BD_rateT, 0)) = array2table(nan);
% BD_snrT(ismember(BD_snrT, 0)) = array2table(nan);

fprintf("BD Rate : %s \n", metric_name);
disp(BD_rateT)
fprintf("BD SNR : %s \n", metric_name);
disp(BD_snrT)
end