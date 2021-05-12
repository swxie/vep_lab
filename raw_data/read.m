%一个用于从该文件夹下读取原始文本数据，转化为mat，并存入data文件夹中的脚本
clear;
clc;
%文件的第一列为激励，第二、三列为脑电，第四列为眼电
%设置参数
name = '.\*_*';
filelist = dir(name);
file_num  = size(filelist,1);
%设置滤波参数
fft_left = 0.8;%带通滤波器的低截止频率
fft_right = 40;%带通滤波器的高截止频率
filt_order = 2;%滤波器阶数
Fs = 1e5;%降采样前的频率
fs  = 100;%降采样后的频率
t_dis = 1;%舍弃的头尾片段，避免开始测量和
%设置片段参数
delta_t = 1;%刺激的时间间隔
window_length = delta_t * fs;%刺激间隔对应的窗口长度
interest_length = 0.3 * fs;%感兴趣区间长度
%读取文件循环
for i = 1:file_num
    %读取文件
    data = readtable(strcat('.\', filelist(i).name));
    data = table2array(data);
    %去除可能的空值
    data = rmmissing(data);
    %滤波:50Hz以及带通滤波
    [a,b] = butter(filt_order, [45, 55]  / (Fs / 2), "stop");
    data(:,2) = filtfilt(a, b, data(:,2));
    data(:,3) = filtfilt(a, b, data(:,3));
    data(:,4) = filtfilt(a, b, data(:,4));
    f_fir = [fft_left, fft_right];
    [a,b] = butter(filt_order, f_fir * 2 / Fs, "bandpass");
    data(:,2) = filtfilt(a, b, data(:,2));
    data(:,3) = filtfilt(a, b, data(:,3));
    data(:,4) = filtfilt(a, b, data(:,4));
    %降采样和参数修正
    data = [resample(data(:,1), fs, Fs), resample(data(:,2), fs, Fs), resample(data(:,3), fs, Fs),resample(data(:,4), fs, Fs)];
    data = data(t_dis * fs + 1 : end - t_dis * fs , :);
    save(strcat('..\data\',filelist(i).name, '.mat'), 'data');
end

