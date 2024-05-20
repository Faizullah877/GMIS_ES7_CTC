function [residue8bit, idx, idx_bytes] = convert_10bit_residue_to_8bit(residue, compressed_set_folder, imgNo)

index_dir = fullfile(compressed_set_folder, "IdxImages");
if(~exist(index_dir, 'dir')), mkdir(index_dir); end
index_images_dir = fullfile(index_dir, sprintf("IdxImages_IMG%03d",imgNo));
if(~exist(index_images_dir, 'dir')), mkdir(index_images_dir); end

tmp = (int16(residue))+128; %now tmp is [-128 -> 382] from [-255 -> +255]
diff_C1 = tmp(:,:,1);
diff_C2 = tmp(:,:,2);
diff_C3 = tmp(:,:,3);
idx = struct();
ss = size(diff_C1);
y1_C1 = zeros(ss);
y1_C2 = zeros(ss);
y1_C3 = zeros(ss);
i1_C1 = zeros(ss);
i1_C2 = zeros(ss);
i1_C3 = zeros(ss);
i2_C1 = zeros(ss);
i2_C2 = zeros(ss);
i2_C3 = zeros(ss);

[h, w, ~] = size(tmp);
for m=1:h
    for n=1:w
        if(diff_C1(m,n)<=255 && diff_C1(m,n)>0)
            y1_C1(m,n)=diff_C1(m,n);
        elseif(diff_C1(m,n)>255)
            y1_C1(m,n) = 510 - diff_C1(m,n);
            i1_C1(m,n) = 1;
        else
            y1_C1(m,n) = -(diff_C1(m,n));
            i2_C1(m,n) = 1;
        end
    end
end
for m=1:h
    for n=1:w
        if(diff_C2(m,n)<=255 && diff_C2(m,n)>0)
            y1_C2(m,n)=diff_C2(m,n);
        elseif(diff_C2(m,n)>255)
            y1_C2(m,n) = 510 - diff_C2(m,n);
            i1_C2(m,n) = 1;
        else
            y1_C2(m,n) = -(diff_C2(m,n));
            i2_C2(m,n) = 1;
        end
    end
end
for m=1:h
    for n=1:w
        if(diff_C3(m,n)<=255 && diff_C3(m,n)>0)
            y1_C3(m,n)=diff_C3(m,n);
        elseif(diff_C3(m,n)>255)
            y1_C3(m,n) = 510 - diff_C3(m,n);
            i1_C3(m,n) = 1;
        else
            y1_C3(m,n) = -(diff_C3(m,n));
            i2_C3(m,n) = 1;
        end
    end
end

residue8bit = cat(3, uint8(y1_C1), uint8(y1_C2), uint8(y1_C3));
imwrite(logical(i1_C1), fullfile(index_images_dir, sprintf("index1_C1_%d.png", imgNo)));
imwrite(logical(i1_C2), fullfile(index_images_dir, sprintf("index1_C2_%d.png", imgNo)));
imwrite(logical(i1_C3), fullfile(index_images_dir, sprintf("index1_C3_%d.png", imgNo)));
imwrite(logical(i2_C1), fullfile(index_images_dir, sprintf("index2_C1_%d.png", imgNo)));
imwrite(logical(i2_C2), fullfile(index_images_dir, sprintf("index2_C2_%d.png", imgNo)));
imwrite(logical(i2_C3), fullfile(index_images_dir, sprintf("index2_C3_%d.png", imgNo)));

idx_info = dir(index_images_dir);
idx_bytes = 0;
for i =1: numel(idx_info)
    idx_bytes = idx_bytes + idx_info(i).bytes;
end


idx.i1_C1 = i1_C1;
idx.i1_C2 = i1_C2;
idx.i1_C3 = i1_C3;
idx.i2_C1 = i2_C1;
idx.i2_C2 = i2_C2;
idx.i2_C3 = i2_C3;
end

