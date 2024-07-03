function [residue8bit, idx1, idx2] = sju_residue_dpcm_pixels_forward(previous, next)
diff = int16(previous) - int16(next);
tmp = diff+128; %now tmp is [-128 -> 382] from [-255 -> +255]
diff_C1 = tmp(:,:,1);
diff_C2 = tmp(:,:,2);
diff_C3 = tmp(:,:,3);
ss = size(diff_C1);
y1_C1 = zeros(ss);
y1_C2 = zeros(ss);
y1_C3 = zeros(ss);
i1_C1 = zeros(ss, 'logical');
i1_C2 = zeros(ss, 'logical');
i1_C3 = zeros(ss, 'logical');
i2_C1 = zeros(ss, 'logical');
i2_C2 = zeros(ss, 'logical');
i2_C3 = zeros(ss, 'logical');

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
idx1 = cat(3, i1_C1, i1_C2, i1_C3);
idx2 = cat(3, i2_C1, i2_C2, i2_C3);
end