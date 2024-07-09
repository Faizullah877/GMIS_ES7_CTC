clc
clear all

orig_res_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_5\tmp\residue_002.ppm";
rec_rec_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_5\tmp\dec_image_002.ppm";
com_res_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_5\Compressed\seq_002.jpg";
res_info = dir(com_res_file);
sizeb = res_info.bytes;

orig_img = imread(orig_res_file);
rec_img = imread(rec_rec_file);
no_rows = size(orig_img,1);
no_cols = size(orig_img,2);
img_ncolors = size(orig_img,3);
img_bpp = (sizeb*8) / (no_cols*no_rows);
psnvra = psnr(rec_img, orig_img);
fprintf("Image BPP : %f => PSNR = %f\n", img_bpp, psnvra);

orig_res_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_90\tmp\residue_002.ppm";
rec_rec_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_90\tmp\dec_image_002.ppm";
com_res_file = "D:\GMIS_EXPs\Evaluate_10Fv2\GolfClubHead_CityScape_10thF\SJU_Arch_DPCM_PIXELSv2\q_90\Compressed\seq_002.jpg";
res_info = dir(com_res_file);
sizeb = res_info.bytes;

orig_img = imread(orig_res_file);
rec_img = imread(rec_rec_file);
no_rows = size(orig_img,1);
no_cols = size(orig_img,2);
img_ncolors = size(orig_img,3);
img_bpp = (sizeb*8) / (no_cols*no_rows);
psnvra = psnr(rec_img, orig_img);
fprintf("Image BPP : %f => PSNR = %f\n", img_bpp, psnvra);

