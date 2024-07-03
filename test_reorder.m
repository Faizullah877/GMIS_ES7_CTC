clc
clear all

arr1 = 0 : 1 : 255;
arr2 = zeros(256);
for i = 1: 64
    for j=1:4
        arr2(i + 64 *j ) = arr1();
    end
end
disp(arr2)