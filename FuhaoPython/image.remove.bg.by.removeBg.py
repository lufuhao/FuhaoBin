#!/usr/bin/env python3
from removebg import RemoveBg
 
rmbg = RemoveBg('XN6zNsMKmZoeZLiTvtyCUdEF', 'error.log')
# 第一个参数是 API，第二个参数是将错误输出到日志文件
rmbg.remove_background_from_img_file('girl.jpg') # 括号内是图片地址

