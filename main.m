%����������ȡ�ļ������м��㲢��ӡ�ź����������������Ա�ͼ�Ľű�
filelist = dir('.\data\*_1.mat');
isPlot = 1;%��ͼ
for i = 1 : length(filelist)
    name = filelist(i).name;
    name(end - 4) = '?';
    print_name = name(1 : end - 6);
    name = sprintf('.\\data\\%s', name);
    %�Զ�������������Գ���ʹ�ò�ͬ��arg�����д���
    if i == 1
        fprintf('����\t����ǰ\t�����\t���\n');
    end
    arg = init_arg();
    arg.dwt_order = 3;
    arg.threshold = 0.1;
    arg.slice = 1;
    arg.wave = "db2";
    [xcorr_after, xcorr_ori] = testInnerPatient(name, isPlot, arg);
    change = xcorr_after - xcorr_ori;
    fprintf("%s\t%.2f\t%.2f\t%.2f\n", print_name, xcorr_ori, xcorr_after, change);
    %�Զ������������
end