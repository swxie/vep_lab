%�������ݼ�����֤��
train_set = ["lle", "st5", "wyx","xsw", "ybf","zxy"];
test_set = ["st1", "st2", "st3", "st4"];
cur_set = train_set; %ѡ��ʹ�õļ���

% %����1 �� ��svm�������ж���������,��ʱ��ʹ��С���任
% arg = init_arg();
% arg.dwt_order = 0;
% dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
% nu = [];
% sigma = [];
% acc = [];
% nu_exp_range = 0 : -0.5 : -3;%nu��ָ����Χ
% sigma_exp_range = -3 : 0.5 : 3;%sigma��ָ����Χ
% for i = nu_exp_range
%     for j = sigma_exp_range
%         nu = [nu; 10^i];
%         sigma = [sigma; 10^j];
%         acc_cur = [];
%         for  k = 1 : length(cur_set)
%             [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1, 10^i, 10^j);
%             acc_cur = [acc_cur; mean_acc];
%         end
%         acc = [acc; mean(acc_cur)];
%     end
% end
% nu = reshape(nu, length(sigma_exp_range), length(nu_exp_range))';
% acc = reshape(acc, length(sigma_exp_range), length(nu_exp_range))';
% sigma = reshape(sigma, length(sigma_exp_range), length(nu_exp_range))';
% figure;
% pcolor(nu, sigma, acc);
% set(gca','yscale','log');
% set(gca','xscale','log');
% xlabel("nu");
% ylabel("sigma");
% shading interp;
% colorbar;
% %����1����


% %����2:��ϸ����������
% arg = init_arg();
% arg.dwt_order = 0;
% dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
% nu = [];
% sigma = [];
% acc = [];
% nu_range = 2e-1 : 2e-1 : 1;
% sigma_range = 1e-2 : 1e-2 : 5e-2;
% for i = nu_range
%     for j = sigma_range
%         nu = [nu; i];
%         sigma = [sigma; j];
%         acc_cur = [];
%         for  k = 1 : length(cur_set)
%             [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1, i, j);
%             acc_cur = [acc_cur; mean_acc];
%         end
%         acc = [acc; mean(acc_cur)];
%     end
% end
% nu = reshape(nu, length(sigma_range), length(nu_range))';
% acc = reshape(acc, length(sigma_range), length(nu_range))';
% sigma = reshape(sigma, length(sigma_range), length(nu_range))';
% figure;
% pcolor(nu, sigma, acc);
% xlabel("nu");
% ylabel("sigma");
% shading interp;
% colorbar;
% %����2����


% %����3�����Բ�ͬС���任�����Խ����Ӱ��
% arg = init_arg();
% arg.dwt_order = 0;
% db = [];
% order = [];
% result = [];
% %%�����ޱ任���
% db = [db; ""];
% order = [order; 0];
% cur_result = [];
% for  k = 1 : length(cur_set)
%     dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
%     [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1);
%     cur_result = [cur_result, mean_acc];
% end
% cur_result = [cur_result, mean(cur_result)];
% result = [result; cur_result];
% %%i��j���Ը������ѡ��
% for i = 1 : 3
%     for j = 1 : 3
%         arg.dwt_order = i;
%         arg.wave = "db" + j;
%         db = [db; arg.wave];
%         order = [order; arg.dwt_order];
%         cur_result = [];
%         for  k = 1 : length(cur_set)
%             dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
%             [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1, 0.5, 0.04);
%             cur_result = [cur_result, mean_acc];            
%         end
%         cur_result = [cur_result, mean(cur_result)];
%         result = [result; cur_result];
%     end
% end
% fprintf("С����������");
% fprintf("����\tС��\t׼ȷ��\n");
% for i = 1 : size(result, 1)
%    fprintf("%d\t%s\t%.3f\n", order(i), db(i), result(i, end));
% end
% %����3����


% %����4�� ������Լ����
% cur_set = test_set;
% cur_result = [];
% arg = init_arg();
% arg.dwt_order = 3;
% arg.wave = "db2";
% for  k = 1 : length(cur_set)
%     dataSet = make_dataSet(".\data\" + cur_set(k) + "*", arg);
%     [mean_acc, std_acc, entropy] = svm_classify(dataSet, 10, 0.1, 0.5, 0.04);
%     cur_result = [cur_result, mean_acc];
% end
% display(cur_result);
% %����4����

%�����������ݼ��ĺ���
function [dataSet] = make_dataSet(name, arg)
%�����������ݼ��ĺ�����nameΪ����mat�ļ���·����argΪ�����ṹ�壬��ѡ
if nargin == 1
    arg = init_arg();
end
%�����ļ�����
filelist = dir(name);
file_num  = size(filelist,1);
fs  = arg.fs;
%����Ƭ�β���
delta_t = arg.delta_t;
window_length = delta_t * fs;
interest_length = arg.interest_length * fs;
%����С���任����
dwt_order = arg.dwt_order;
wave = arg.wave;
%�������ݼ�����
data_set_pos = [];
data_set_neg = [];
%��ʷ��������
sum_count = 1;
%��ȡ�ļ�ѭ��
for i = 1:file_num
    %��ȡ�ļ�
    load(strcat('./data/', filelist(i).name), 'data');
    %���켤������
    buffer = data(1:window_length, 1);
    diff = conv(buffer,[1;-1]);
    diff(1:5) = zeros(1,5);
    diff(end - 4 : end) = zeros(1,5);
    max_index = find(diff == max(diff),1);
    end_mat = (max_index + interest_length - 1) : window_length : size(data,1);
    begin_mat = end_mat - interest_length + 1;
    %����ƽ��
    count = floor(length(begin_mat) / sum_count);
    data = reshape(data(1 : window_length * count * sum_count, :), sum_count, count * window_length, size(data, 2) );
    count = count - 1;
    if (sum_count ~= 1)
        data = mean(data);
    end
    data = reshape(data, size(data, 2), size(data, 3));
    begin_mat = round(begin_mat);
    end_mat = round(end_mat);
    %����vep����
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
%���ؽ��
dataSet = struct("data_set_pos", data_set_pos, "data_set_neg", data_set_neg, ...
"label_neg", label_neg, "label_pos", label_pos);
%save("dataSet.mat","data_set_pos","data_set_neg","label_neg","label_pos");
%state = 1;
end

function [mean_acc, std_acc, mean_entropy] = svm_classify(dataSet, k, pos_ratio, nu, sigma)
%����k��svm���࣬�����������䵽ԭ����pos_ratio��
if nargin < 4
    nu = 0.5;
    sigma = 1;
end
%��ȡ����
data_set_pos = dataSet.data_set_pos;
data_set_neg = dataSet.data_set_neg;
label_pos = dataSet.label_pos;
label_neg = dataSet.label_neg;
%load('dataSet.mat');
%���ɽ������
acc_test = zeros(k, 1);
entropy_test = zeros(k, 1);
%��ʼѭ��
for i = 1 : k
    %�������ݼ�
    data_set_neg = data_set_neg(randperm(size(data_set_neg, 1)), :);%����˳��
    count_neg = round(size(data_set_pos, 1) * pos_ratio);
    train_times = floor(size(data_set_neg, 1) / count_neg);
    %���ɽ������
    acc = zeros(train_times, 1);
    entropy = zeros(train_times, 1);
    %��ʼѭ��
    for j = 1 : train_times
        %����ѵ����
        begin_index = (j - 1) * count_neg + 1;
        end_index = j * count_neg;
        data_set = [data_set_pos; data_set_neg(begin_index : end_index,:)];
        true_label = [label_pos; label_neg(begin_index : end_index)];
        label = ones(size(data_set, 1),1);
        %ѵ��һ����SVM��
        SVMModel = fitcsvm(data_set, label,'kernelScale',sigma, ...
            'nu', nu);
        %������ȷ��
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


