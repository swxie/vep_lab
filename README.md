# vep_lab
一个用于实践、测试FVEP信号增强算法对FVEP信号可重复性指标影响的代码框架
## 文件组成
- raw_data目录：原始字符文件数据，且包含了一个用于从字符文件中读取数据并存为data目录下的二进制文件的脚本read.m
- data目录：存放读取后的二进制文件（.mat）数据，内容为6e6×4的60s四路FVEP原始数据，文件名皆为[受试者]_i的形式，其中i必须从1开始，且必须小于10
- **main.m: 主运行脚本， 用于从data目录中批量读取mat，进行信号处理计算并打印处理前后信号质量分数的变化，并画出对比图**
- **classifier_test.m：一个研究信号处理算法的示例，包含了四个用来研究分类器性能、进行调优的函数，可以通过去注释来选择调用哪一段代码。**
- signal_enhance.m：信号处理函数的主文件，包括了滤波和降采样、小波变换、vep分段、去除低质量片段等步骤
- init_arg.m: 包含了用于初始化信号处理参数结构体arg的函数
- eegICA.m：ICA相关函数，没有被main.m使用，可以自选加入
- testInnerPatient.m：测试信号处理性能的函数，输入为测试信号文件的通配符表达式，返回处理后和处理前信号的可重复性指标（互相关），并画出处理前后的对比图
- README.md：本文件
## 任务
- 将源数据放入对应的目录
  - 形如"xxx_x.mat"的二进制mat文件应放入data目录下
  - 形如"xxx_x"的文本文件，应放入raw_data目录下，并用read.m脚本产生对应的data目录文件
- 阅读signal_enhance代码，了解整个系统的数据流
- 修改signal_enhance程序，来尝试不同的信号处理算法/信号处理参数，然后运行main.m脚本来观察不同信号处理算法/信号处理参数的效果差别
- 你可以在main.m脚本中，通过修改传入的arg结构体来修改信号处理参数，进而实现对参数的自动化测试
- 你也可以编写自己的脚本，基于数据做一些自由的研究，classifier_test.m会是一个合适的例子

## 数据链接
南大学生可以打开NJU云盘链接：https://box.nju.edu.cn/d/cbc6762a0e054890b1aa/
