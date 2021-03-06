%用于批量读取文件，进行计算并打印信号质量分数、画出对比图的脚本
filelist = dir('.\data\*_1.mat');
isPlot = 1;%画图
for i = 1 : length(filelist)
    name = filelist(i).name;
    name(end - 4) = '?';
    print_name = name(1 : end - 6);
    name = sprintf('.\\data\\%s', name);
    %自定义代码区，可以尝试使用不同的arg来进行处理
    if i == 1
        fprintf('样本\t处理前\t处理后\t提高\n');
    end
    arg = init_arg();
    arg.dwt_order = 3;
    arg.threshold = 0.1;
    arg.slice = 1;
    arg.wave = "db2";
    [xcorr_after, xcorr_ori] = testInnerPatient(name, isPlot, arg);
    change = xcorr_after - xcorr_ori;
    fprintf("%s\t%.2f\t%.2f\t%.2f\n", print_name, xcorr_ori, xcorr_after, change);
    %自定义代码区结束
end