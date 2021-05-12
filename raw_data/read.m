%һ�����ڴӸ��ļ����¶�ȡԭʼ�ı����ݣ�ת��Ϊmat��������data�ļ����еĽű�
clear;
clc;
%�ļ��ĵ�һ��Ϊ�������ڶ�������Ϊ�Ե磬������Ϊ�۵�
%���ò���
name = '.\*_*';
filelist = dir(name);
file_num  = size(filelist,1);
%�����˲�����
fft_left = 0.8;%��ͨ�˲����ĵͽ�ֹƵ��
fft_right = 40;%��ͨ�˲����ĸ߽�ֹƵ��
filt_order = 2;%�˲�������
Fs = 1e5;%������ǰ��Ƶ��
fs  = 100;%���������Ƶ��
t_dis = 1;%������ͷβƬ�Σ����⿪ʼ������
%����Ƭ�β���
delta_t = 1;%�̼���ʱ����
window_length = delta_t * fs;%�̼������Ӧ�Ĵ��ڳ���
interest_length = 0.3 * fs;%����Ȥ���䳤��
%��ȡ�ļ�ѭ��
for i = 1:file_num
    %��ȡ�ļ�
    data = readtable(strcat('.\', filelist(i).name));
    data = table2array(data);
    %ȥ�����ܵĿ�ֵ
    data = rmmissing(data);
    %�˲�:50Hz�Լ���ͨ�˲�
    [a,b] = butter(filt_order, [45, 55]  / (Fs / 2), "stop");
    data(:,2) = filtfilt(a, b, data(:,2));
    data(:,3) = filtfilt(a, b, data(:,3));
    data(:,4) = filtfilt(a, b, data(:,4));
    f_fir = [fft_left, fft_right];
    [a,b] = butter(filt_order, f_fir * 2 / Fs, "bandpass");
    data(:,2) = filtfilt(a, b, data(:,2));
    data(:,3) = filtfilt(a, b, data(:,3));
    data(:,4) = filtfilt(a, b, data(:,4));
    %�������Ͳ�������
    data = [resample(data(:,1), fs, Fs), resample(data(:,2), fs, Fs), resample(data(:,3), fs, Fs),resample(data(:,4), fs, Fs)];
    data = data(t_dis * fs + 1 : end - t_dis * fs , :);
    save(strcat('..\data\',filelist(i).name, '.mat'), 'data');
end

