function [left_vep, right_vep, learning_curve, left_vep_before, right_vep_before] = signal_enhance(data, arg)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
%����vep��·�ź����룬���غϳɵ�left_vep��right_vep
%   ����ȫ�ֱ���
if nargin == 1
    arg = init_arg();
end
%�����˲�����
fft_left = arg.fft_left;%��ͨ�˲����ĵͽ�ֹƵ��
fft_right = arg.fft_right;%��ͨ�˲����ĸ߽�ֹƵ��
filt_order = arg.filt_order;%�˲�������
Fs = arg.Fs;%������ǰ��Ƶ��
fs  = arg.fs;%���������Ƶ��
t_dis = arg.t_dis;%������ͷβƬ�Σ����⿪ʼ������
%����Ƭ�β���
delta_t = arg.delta_t;%�̼���ʱ����
window_length = delta_t * fs;%�̼������Ӧ�Ĵ��ڳ���
interest_length = arg.interest_length * fs;%����Ȥ���䳤��
%����С���任����
dwt_order = arg.dwt_order;%dwt����
wave = arg.wave;%С������
%ɸѡ��ز���
slice = arg.slice;%ÿ��ɸѡɸѡ���ĵ������ź�Ƭ�θ���
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
%����ֶ�
buffer = data(1:window_length, 1);
diff = conv(buffer,[1;-1]);
diff(1:5) = zeros(1,5);
diff(end - 4 : end) = zeros(1,5);
max_index = find(diff == max(diff),1);
end_mat = (max_index + interest_length - 1) : window_length : size(data,1);
begin_mat = end_mat - interest_length + 1;
left_vep = zeros(interest_length, length(begin_mat));
right_vep = zeros(interest_length, length(begin_mat));
%�����������
for i = 1:length(begin_mat)
    left_vep(:, i) = data(begin_mat(i) : end_mat(i), 2);
    right_vep(:, i) = data(begin_mat(i) : end_mat(i), 3);
end
left_vep_before = left_vep;
right_vep_before = right_vep;
%%����dwt�󳤶�
temp = left_vep(:,1);
for k = 1 : dwt_order
    temp = dwt(temp,wave);
end
dwt_length = length(temp);
%%��ʼ��vep����
vep = zeros(2 * dwt_length, length(begin_mat));
for i = 1:length(begin_mat)
    left_temp = left_vep(:, i);
    right_temp = right_vep(:, i);
    for k = 1 : dwt_order
        left_temp = dwt(left_temp,wave);
        right_temp = dwt(right_temp, wave);
    end
    vep(:, i) = [left_temp; right_temp];
end
% figure;
% subplot(2, 2, 1);
% plot(sum(left_vep, 2));
% subplot(2, 2, 2);
% plot(sum(right_vep, 2));
%%���е����б�ɾ��
learning_curve = (mean(mean(corrcoef(left_vep))) + mean(mean(corrcoef(right_vep)))) / 2;
while 1
    if mod(size(left_vep, 2), slice) == 0 
        mask = svm_classify(vep, slice);
    else
        mask = svm_classify(vep, mod(size(left_vep, 2), slice));
    end
    if sum(mask) == length(mask)
        break;
    end
    left_vep_new = left_vep(:, mask > 0);
    right_vep_new = right_vep(:, mask > 0);
    left_vep_cur = zeros(size(left_vep_new, 1), slice);
    right_vep_cur = zeros(size(right_vep_new, 1), slice);
    part_length = round(size(left_vep_new, 2) / slice);
    for i = 1 : slice
        left_vep_cur(:, i) = mean(left_vep(:, (i - 1) * part_length + 1 : i * part_length), 2);
        right_vep_cur(:, i) = mean(right_vep(:, (i - 1) * part_length + 1 : i * part_length), 2);
    end
    learning_curve = [learning_curve; (mean(mean(corrcoef(left_vep_cur))) + mean(mean(corrcoef(right_vep_cur)))) / 2];
    if learning_curve(end) < learning_curve(end - 1)
        break;
    end
    left_vep = left_vep_new;
    right_vep = right_vep_new;
    vep = vep(:, mask > 0);
end
% subplot(2, 2, 3);
% plot(sum(left_vep, 2));
% subplot(2, 2, 4);
% plot(sum(right_vep, 2));
end

function [mask] = svm_classify(data, slice)
    SVMModel = fitcsvm(data', ones(size(data, 2), 1),'kernelScale','auto', ...
            'OutlierFraction',0.05);
    [~,scorePred] = predict(SVMModel, data');
    t = sort(scorePred);
    if length(scorePred) < slice + 1
        mask = ones(size(scorePred));
    else
        mask = scorePred >= t(slice + 1);
    end
end

