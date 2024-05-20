% convert from rgb to xyb using reference C code
function [outxyb] = conv_xyb_to_srgb(inxyb)
    addpath('pfm');

    pfm_folder = 'data\';
    prefix_name = 'inXYB';
    in_pfm_ref_name = [pfm_folder prefix_name '.pfm'];
    write_pfm(inxyb, in_pfm_ref_name, 'l', 0);
    
    
    % read pfm and convert to RGB and save as ppm
    prefix_name = 'outRGB';
    out_ppm_ref_name = [pfm_folder prefix_name '.ppm'];
    
    EXE_DIR = 'C:\Users\reji_\MatlabProjects\xyb_colorspace\xyb_converter-main\xybConversion\x64\Debug\';
    CONV_XYB = ['rgbConversion.exe ' in_pfm_ref_name  ' ' out_ppm_ref_name];
    EXE_COMMAND = [EXE_DIR CONV_XYB];
    [status_conv,cmdout_conv] = system(EXE_COMMAND);

    outxyb = double(imread(out_ppm_ref_name));
end