#!/usr/bin/env python3

#https://www.paddlepaddle.org.cn/



################ 1. install paddlepaddle ######################
#python3 -m pip install paddlepaddle -i https://mirror.baidu.com/pypi/simple
###needs pathlib, click, joblib, regex, tqdm, nltk, gast, rarfile, pyyaml, funcsigs, paddlepaddle

###python环境测试是否安装成功
import paddle.fluid
paddle.fluid.install_check.run_check()
###如果能看到 Your Paddle Fluid is installed successfully 就表示安装成功了。



################ 2. install paddlehub ######################
###要实现本文的一键批量扣图需求，需要借助PaddleHub人像分割模型来实现。
###https://www.paddlepaddle.org.cn/hub
###https://github.com/PaddlePaddle/PaddleHub
###更多PaddleHub预训练模型教程合集课程可见：
###https://aistudio.baidu.com/aistudio/course/introduce/1070
###介绍完了项目，接下来我们开始在线安装 paddlehub ：
#pip install -i https://mirror.baidu.com/pypi/simple paddlehub
###或者按指定版本安装：
#pip install paddlehub==1.6.0 -i https://pypi.tuna.tsinghua.edu.cn/simple



################ 3. 一键扣图代码实现 ######################
###我们的实现步骤很简单：
###    导入模块
###    加载模型
###    获取图片文件
###    调用模块抠图
###其中扣图功能主要采用PaddleHub DeepLabv3+模型deeplabv3p_xception65_humanseg。
###下面我们看具体扣图代码实现(demo.py)：

import os
import paddlehub as hub

# 加载模型
humanseg = hub.Module(name='deeplabv3p_xception65_humanseg')  
base_dir = os.path.abspath(os.path.dirname(__file__))

# 获取当前文件目录
path = os.path.join(base_dir, 'images/')
# 获取文件列表
files = [path + i for i in os.listdir(path)]  
print(files)
# 抠图
results = humanseg.segmentation(data={'image': files})  
for result in results:
    print(result)
###示例中，我将图片放在代码文件夹的同级目录 images文件夹下，运行代码后，输出的抠图图片会自动放在代码同级目录的 humanseg_output 目录下，文件名称跟原图片的名称相同，但是文件格式是 png 。



################ 4. 需要注意的坑 ################
###在运行示例代码时，如果没有单独安装模型deeplabv3p_xception65_humanseg，默认会自动在执行前进行安装。但安装完成后，执行结果并没有生成扣图结果及humanseg_output目录，输出结果类似如下所示：
###正常情况下，在生成扣图数据，打印results时，应该是类似如下结构才对：
###可以通过单独安装模型并指定安装版本来解决。
#hub install deeplabv3p_xception65_humanseg==1.0.0
###具体原因没有细究，默认自动安装模型时，版本为1.2.0，猜测由于还是模型版本不兼容问题导致。


################ 5. 总结 ################
###本文基于 paddlepaddle 平台，利用PaddleHub DeepLabv3+模型(deeplabv3p_xception65_humanseg)，使用简单的五行代码就实现了批量抠图。有些读者可能会想，上述示例中提供的代码行数不止五行代码吧，在上述示例中，真正实现扣图的主代码其实只需要下面五行：
#humanseg = hub.Module(name='deeplabv3p_xception65_humanseg')  
#base_dir = os.path.abspath(os.path.dirname(__file__))
#path = os.path.join(base_dir, 'images/')
#files = [path + i for i in os.listdir(path)]  
#results = humanseg.segmentation(data={'image': files})  
###利用PaddleHub DeepLabv3+模型 不仅可以实现一键扣图，还可以进行图片合成，视频合成等。利用好它不仅解放了人的双手和双眼，而且为某些程序猿/程序媛的装逼工具箱提供了一件宝器。下次如果碰到某个女生或者闺蜜在为抠图发愁，别忘了掏出神器，赢得芳心哦！
###paddlepaddle作为一款开源的深度学习平台，本文介绍的扣图训练模型只是其中的冰山一角，实战训练预测模型种类还远远不止，更多的场景结合，读者们可自行挖掘。
