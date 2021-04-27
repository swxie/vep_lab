%一个用于从该文件夹下读取原始文本数据，转化为mat，并存入data文件夹中的脚本
clear;
clc;
%文件的第一列为激励，第二列为脑电
%设置参数
name = '.\*_*';
filelist = dir(name);
file_num  = size(filelist,1);
%读取文件循环
for i = 1:file_num
    %读取文件
    data = readtable(strcat('.\', filelist(i).name));
    data = table2array(data);
    save(strcat('..\data\',filelist(i).name, '.mat'), 'data');
end

