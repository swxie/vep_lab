%区分数据集和验证集
train_set = ["liule", "stu5", "wyx","xsw", "ybf","zxy"];
test_set = ["stu1", "stu2", "stu3", "stu4"];
cur_set = train_set; %选择使用的集合
%测试不同小波变换参数对结果的影响
fprintf("无小波变换:\n");
fprintf("序号\t平均准确率\t标准差\t熵\n");
arg = init_arg();
arg.dwt_order = 0;
for  k = 1 : length(cur_set)
    dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
    [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1);
    fprintf("%s\t%f\t%f\t%f\n", cur_set(k), mean_acc, std_acc,entropy);
end
for i = 1 : 3
    for j = 1 : 3
        arg.dwt_order = j;
        arg.wave = "db" + i;
        fprintf("%s小波和%d次变换:\n", arg.wave, arg.dwt_order);
        fprintf("序号\t平均准确率\t标准差\t熵\n");
        for  k = 1 : length(cur_set)
            dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
            [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1);
            fprintf("%s\t%f\t%f\t%f\n", cur_set(k), mean_acc, std_acc, entropy);
        end
    end
end
%对svm参数进行对数粗搜索
arg = init_arg();
dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
nu = [];
sigma = [];
acc = [];
nu_exp_range = 0 : -0.5 : -3;%nu的指数范围
sigma_exp_range = -3 : 0.5 : 3;%sigma的指数范围
for i = nu_exp_range
    for j = sigma_exp_range
        nu = [nu; 10^i];
        sigma = [sigma; 10^j];
        acc_cur = [];
        for  k = 1 : length(cur_set)
            [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1, 10^i, 10^j);
            acc_cur = [acc_cur; mean_acc];
        end
        acc = [acc; mean(acc_cur)];
    end
end
nu = reshape(nu, length(sigma_exp_range), length(nu_exp_range))';
acc = reshape(acc, length(sigma_exp_range), length(nu_exp_range))';
sigma = reshape(sigma, length(sigma_exp_range), length(nu_exp_range))';
figure;
pcolor(nu, sigma, acc);
set(gca','yscale','log');
set(gca','xscale','log');
xlabel("nu");
ylabel("sigma");
shading interp;
colorbar;


%用于生成数据集的函数
function [dataSet] = make_dataSet(name, arg)
%用于生成数据集的函数，name为数据mat文件的路径，arg为参数结构体，可选
if nargin == 1
    arg = init_arg();
end
%设置文件参数
filelist = dir(name);
file_num  = size(filelist,1);
Fs = arg.Fs;
%设置滤波参数
fft_left = arg.fft_left;
fft_right = arg.fft_right;
filt_order = arg.filt_order;
fs  = arg.fs;
t_dis = arg.t_dis;
%设置片段参数
delta_t = arg.delta_t;
window_length = delta_t * fs;
interest_length = arg.interest_length * fs;
%设置小波变换参数
dwt_order = arg.dwt_order;
wave = arg.wave;
%创建数据集容器
data_set_pos = [];
data_set_neg = [];
%历史遗留参数
sum_count = 1;
%读取文件循环
for i = 1:file_num
    %读取文件
    load(strcat('./data/', filelist(i).name), 'data');
    %滤波
    f_fir = [fft_left, fft_right];
    [a,b] = butter(filt_order, f_fir * 2 / Fs, "bandpass");
    data(:,2) = filtfilt(a, b, data(:,2));
    data(:,3) = filtfilt(a, b, data(:,3));
    data(:,4) = filtfilt(a, b, data(:,4));
    %降采样和参数修正
    data = [resample(data(:,1), fs, Fs), resample(data(:,2), fs, Fs), resample(data(:,3), fs, Fs),resample(data(:,4), fs, Fs)];
    data = data(t_dis * fs + 1 : end - t_dis * fs , :);
    %构造激励矩阵
    buffer = data(1:window_length, 1);
    diff = conv(buffer,[1;-1]);
    diff(1:5) = zeros(1,5);
    diff(end - 4 : end) = zeros(1,5);
    max_index = find(diff == max(diff),1);
    end_mat = (max_index + interest_length - 1) : window_length : size(data,1);
    begin_mat = end_mat - interest_length + 1;
    %滑动平均
    count = floor(length(begin_mat) / sum_count);
    data = reshape(data(1 : window_length * count * sum_count, :), sum_count, count * window_length, size(data, 2) );
    count = count - 1;
    if (sum_count ~= 1)
        data = mean(data);
    end
    data = reshape(data, size(data, 2), size(data, 3));
    begin_mat = round(begin_mat);
    end_mat = round(end_mat);
    %构造vep矩阵
    vep = [];
    for j = 1 : count
        %         vep_left_temp = data(begin_mat(j):end_mat(j), 2)';
        %         vep_right_temp = data(begin_mat(j):end_mat(j), 3)';
        %         if (max(vep_left_temp) - min(vep_left_temp) > threshold) || (max(vep_right_temp) - min(vep_right_temp) > threshold)
        %             continue;
        %         end
        data_left = data(begin_mat(j):end_mat(j), 2)';
        data_right = data(begin_mat(j):end_mat(j), 3)';
        data_eye = data(begin_mat(j):end_mat(j), 4)';
        for k = 1 : dwt_order
            data_left = dwt(data_left,wave);
            data_right = dwt(data_right,wave);
            data_eye = dwt(data_eye,wave);
        end
        vep = [vep; data_left, data_right, data_eye];
    end
    if (filelist(i).name(end - 4) == 'b')
        data_set_neg = [data_set_neg;vep];
    else
        data_set_pos = [data_set_pos;vep];
    end
end
data_set_pos = normalize(data_set_pos, 2);
data_set_neg = normalize(data_set_neg, 2);
label_pos = ones(size(data_set_pos, 1),1);
label_neg = zeros(size(data_set_neg, 1),1);
%返回结果
dataSet = struct("data_set_pos", data_set_pos, "data_set_neg", data_set_neg, ...
"label_neg", label_neg, "label_pos", label_pos);
%save("dataSet.mat","data_set_pos","data_set_neg","label_neg","label_pos");
%state = 1;
end

function [mean_acc, std_acc, mean_entropy] = svm_classify(dataSet, k, pos_ratio, nu, sigma)
%进行k次svm分类，正样本数扩充到原来的pos_ratio倍
if nargin < 4
    nu = 0.5;
    sigma = 1;
end
%读取数据
data_set_pos = dataSet.data_set_pos;
data_set_neg = dataSet.data_set_neg;
label_pos = dataSet.label_pos;
label_neg = dataSet.label_neg;
%load('dataSet.mat');
%生成结果数组
acc_test = zeros(k, 1);
entropy_test = zeros(k, 1);
%开始循环
for i = 1 : k
    %生成数据集
    data_set_neg = data_set_neg(randperm(size(data_set_neg, 1)), :);%打乱顺序
    count_neg = round(size(data_set_pos, 1) * pos_ratio);
    train_times = floor(size(data_set_neg, 1) / count_neg);
    %生成结果数组
    acc = zeros(train_times, 1);
    entropy = zeros(train_times, 1);
    %开始循环
    for j = 1 : train_times
        %生成训练集
        begin_index = (j - 1) * count_neg + 1;
        end_index = j * count_neg;
        data_set = [data_set_pos; data_set_neg(begin_index : end_index,:)];
        true_label = [label_pos; label_neg(begin_index : end_index)];
        label = ones(size(data_set, 1),1);
        %训练一分类SVM器
        SVMModel = fitcsvm(data_set, label,'kernelScale',sigma, ...
            'nu', nu);
        %计算正确率
        [~,scorePred] = predict(SVMModel, data_set);
        split_index = round(pos_ratio * length(true_label));
        temp = sort(scorePred);
        threshold = temp(split_index);
        correct_count = size(find(true_label == 0 & (scorePred < threshold) == true_label), 1);
        all_count = size(find(true_label == 0), 1);
        acc(j) = correct_count / all_count;
        scorePred = scorePred - min(scorePred);
        entropy(j) = - mean(log(1 - scorePred(true_label == 0)));
    end
    acc_test(i) = mean(acc);
    entropy_test(i) = mean(entropy);
end
mean_acc = mean(acc_test);
std_acc = std(acc_test);
mean_entropy = mean(entropy_test);
end


