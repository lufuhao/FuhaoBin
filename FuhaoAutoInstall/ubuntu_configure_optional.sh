#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi
abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"



################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20150603

Requirements:
	

Descriptions:
	xxx

Options:
  -h    Print this help message
  -i    CONFIG file
  -t    Number of threads, default: 1
  -s    Not run simulation
  -a    Not run assembly

Example:
  $0 -i ./chr1.fa -t 10

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) opt_s=1;shift 1;;
    -a) opt_a=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0
  else
#    echo "I require $1 but it's not installed.  Aborting." >&2
    echo 1
  fi
#  local cmd=$1
#  if command -v $cmd >/dev/null 2>&1;then
#    echo >&2 $cmd "  :  "`command -v $cmd`
#    exit 0
#  else
#    echo >&2 "Error: require $cmd but it's not installed.  Exiting..."
#    exit 1
#  fi
}
###Usage: array=(`split delimiter string`)
split () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}
#str='ni,hai,a'
#a=(`SplitString ',' $str`)
#echo ${#a[@]} ${a[0]} ${a[1]} ${a[2]}
#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}

### subfunction
function InstallPackages () {
	local -a packagenames=$@
	echo -e "\n\n\n$ProgName is installing $@"
	sudo apt-get install "$@" >> install.log 2>> install.err
}



#################### Command test ###################################
if [ $(CmdExists 'santools') -eq 0 ]; then
	exit 0
else
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'santools') -eq 1 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################




#################### Main ###########################################


### ubuntu firewall
#开启了防火墙，并在系统启 动时自动开启。关闭所有外部对本机的访问，但本机访问外部正常。
sudo ufw enable
sudo ufw default deny
#要查看防火墙的状态
sudo ufw status
#开启／禁用防火墙的某个服务
#sudo ufw allow|deny [service]
#打开或关闭某个端口
#sudo ufw allow smtp　允许所有的外部IP访问本机的25/tcp (smtp)端口
#sudo ufw allow 22/tcp 允许所有的外部IP访问本机的22/tcp (ssh)端口
#sudo ufw allow 53 允许外部访问53端口(tcp/udp)
#sudo ufw allow from 192.168.1.100 允许此IP访问所有的本机端口
#sudo ufw allow proto udp 192.168.0.1 port 53 to 192.168.0.2 port 53
#sudo ufw deny smtp 禁止外部访问smtp服务
#sudo ufw delete allow smtp 删除上面建立的某条规则



#Audacious乱码问题
#播放列表中的文件名显示有些是乱码的, 解决办法如下（经测试有效）
#在audacious上右键选择“首选项”，在“播放列表”中，把标题格式改为“Custom”,再把自定格式改为“%f”（不要引号）



#文件管理器
#类似于Windows下的Total Commander的文件管理器: Krusader，muCommmander，GNOME Commander
InstallPackages krusader
#你可以把Krusader的F2键调用的程序改为Gnome Terminal：
echo "在 设置→配置Krusader→常规 中，“终端”文本框里的“konsole --workdir %d”改为：
gnome-terminal --working-directory %d"



#启用Chrome浏览器的Flash插件
echo "在Chrome浏览器的地址栏里输入：chrome://plugins/
回车进入插件配置界面，你会看到一个名为“libpepflashplayer”和一个名为“Adobe Flash Player”的插件
把前一个插件下方的“始终允许”的复选框勾上，然后点击页面右上方的“详细信息”
在展开的详细信息中，找到“Adobe Flash Player”插件的详细信息处，会发现有两个子项，都叫“Shockwave Flash”，并且两个都是启用的状态。其中有一个的“位置”是以 /usr/lib 开头的，另一个不是，现在不要动那个以 /usr/lib 开头的子项，把另一个“停用”掉
然后重启浏览器，就会发现可以播放页面中的Flash视频了。"



#Ubuntu下好用的计算器软件
sudo apt-get install speedcrunch
sudo apt-get install qalculate


#磁盘管理工具
sudo apt-get install gparted


#让Ubuntu的文件浏览器显示地址栏, 直接按Ctrl+L即可显示


#安装强大的音频处理软件Audacity
sudo apt-get install audacity


#终端（terminal）下，命令前不显示完整路径，而只是显示当前所在的目录名
#编辑当前用户根目录下的 .bashrc 文件，找到以下代码：
#case "$TERM" in
#xterm*|rxvt*)
#	#PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
#	PS1="[\u@ \W]\\$ "
#	;;
#*)
#	;;
#esac


### Teamviewer
sudo apt-get install teamviewer



#当任务栏里的IBus图标丢失时，你可以通过重启IBus来找回：
sudo killall ibus-daemon
sudo ibus-daemon -d
ibus-setup



#在 快速启动列表／快速启动器／Dash主页 中添加一个自定义程序法
cd /usr/share/applications
sudo vim Sublime.desktop

#Version=1.0
#Name=Sublime Text
#Exec=/usr/bin/sublime
#Terminal=false
#Icon=/opt/sublime/Icon/128x128/sublime_text.png
#Type=Application
#Categories=Development



#launcher项的保存位置
#Ubuntu的启动项（desktop entry）是以文件的形式存储在磁盘上的，具体位置就是：
#仅对当前用户有效的： ~/.local/share/applications
#对所有用户都有效的： /usr/share/applications/







#删除自己添加过的软件源: 进入目录 /etc/apt/sources.list.d/，把对应的文件删掉就可以了。



#将旧机器上的Thunderbird配置及数据恢复到新机器上: 
#将旧机器上的 ~/.thunderbird 目录下的“profiles.ini”文件以及一个名字很怪的目录（例如“y88yf66x.default”）拷贝到你新机器上的相同目录下，重启Thunderbird即可。



#Ubuntu下好用的视频编辑软件
#Linux下的视频编辑软件，例如Avidemux，OpenShot，等等。
sudo apt-get install avidemux
sudo apt-get install openshot libavformat-extra-53
#其中，libavformat-extra-53是OpenShot软件在将视频保存为MP4(h.264)格式时依赖的软件包，所以要安装（除非你不把视频导出为这种格式）。



#安装ssh server，使得别的服务器可以ssh连接本机
sudo apt-get install openssh-server
sudo /etc/init.d/ssh restart



#Ubuntu下自带光盘刻录软件: Brasero



###rsync同步本机上的两个目录
#备份文件的话，rsync最方便了：	
rsync -vzrtopgu -progress src_dir/ dest_dir/
#其中，src_dir 是源目录，dest_dir是目标目录。两个目录后都要带/，这样会把目录里的所有文件及子目录都同步到目标目录下。



#apt-get install安装BDB/How to apt-get install BDB(Berkeley DB)
apt-cache search libdb | grep Berkeley
#从列表中找到合适的版本并安装，例如：	
sudo apt-get install libdb5.1 libdb5.1-dev



#配置Ubuntu开机自动连接蓝牙鼠标



#Ubuntu下很好用的视频格式转换软件
sudo apt-get install handbrake



#ubuntu下的epub格式（电子书）阅读软件
sudo apt-get install calibre calibre-bin



#Vim语法高亮
#如果安装好Vim后，发现它对代码文件没有语法高亮（例如，shell脚本）：
#编辑 ~/.vimrc 文件，在里面填上下面的内容：
#if &t_Co > 1
#   syntax enable
#endif



#Ubuntu上的摄像头软件
#已经预装好了，叫“cheese”，中文名是“茄子大头贴”。如果你的系统上没有，就用如下命令安装：
sudo apt-get install cheese


#Ubuntu下好用的看图软件
sudo apt-get install gwenview



#为Ubuntu安装Numix主题和图标
sudo add-apt-repository ppa:numix/ppa
sudo apt-get update
sudo apt-get install numix-gtk-theme numix-icon-theme-circle



#JAVA
#export JAVA_HOME=$RunDir/java/v$jdk_version/$jdk_MachType
#export JRE_HOME=$JAVA_HOME/jre
#export PATH=$JAVA_HOME/bin:$JRE_HOME/bin:$PATH
#export CLASSPATH=.:$JAVA_HOME/lib:$JRE_HOME/lib:$CLASSPATH



if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	exit 1
}

