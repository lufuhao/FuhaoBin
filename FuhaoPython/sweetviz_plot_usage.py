#!/usr/bin/env python3


### https://github.com/fbdesignpro/sweetviz
### sweetviz支持Python 3.6+和Pandas0.25.3+环境
### pip install sweetviz
### sweetviz需要用到基础「os」模块

import sweetviz as sv

### sweetviz使用的原理是，使用一行代码，生成一个数据报告的对象（其中，my_dataframe是pandas中的DataFrame，一种表格型数据结构）

my_report=sv.analyze(my_dataframe)


### default arguments will generate to "SEWEETVIZ_REPORT.html"
my_report.show_html()


### 主函数1
# analyze(source: Union[pd.DataFrame], Tuple[pd.DataFrame, str]], target_feat: str = None, feat_cfg: FeatureConfig=None, pairwise_analysis: str = 'auto')
#数据分析函数中，有4个参数source，target_feat，feat_cfg和pairwise_analysis需要被设置。
#source：以pandas中的DataFrame数据结构、或是DataFrame中的某一类字符串作为分析对象。
#target_feat：需要被标记为目标对象的字符串。
#feat_cfg：需要被跳过、或是需要被强制转换为某种数据类型的特征。
#pairwise_analysis：相关性和其他类型的数据关联可能需要花费较长时间。如果超过了某个阈值，就需要设置这个参数为on或者off，以判断是否需要分析数据相关性。

### 主函数2
#compare()
#my_report=sv.compare([my_dataframe, "Training Data"], [test_df, "Test Data"], "Survived", feature_config)
#如果想要对两个数据集进行对比分析，就使用这个比较函数。
#例子中的my_dataframe和test_df是两个数据集，分别被命名为训练数据和测试数据。
#除了这个被插入的数据集，剩余的参数与analyze中的一致。

### 主函数3
#compare_intra()
#my_report=sv.compare_intra(my_dataframe, my_dataframe["Sex"] == "male", ["Male", "Female"], feature_config)
#想要对数据集中某个栏目下的参数进行分析，就采用这个函数进行。
#例如，如果需要比较“性别”栏目下的“男性”和“女性”，就可以采用这个函数。
#理解这几种函数的变量后，一行代码就能实现Python数据分析。
