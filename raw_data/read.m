%һ�����ڴӸ��ļ����¶�ȡԭʼ�ı����ݣ�ת��Ϊmat��������data�ļ����еĽű�
clear;
clc;
%�ļ��ĵ�һ��Ϊ�������ڶ���Ϊ�Ե�
%���ò���
name = '.\*_*';
filelist = dir(name);
file_num  = size(filelist,1);
%��ȡ�ļ�ѭ��
for i = 1:file_num
    %��ȡ�ļ�
    data = readtable(strcat('.\', filelist(i).name));
    data = table2array(data);
    save(strcat('..\data\',filelist(i).name, '.mat'), 'data');
end

