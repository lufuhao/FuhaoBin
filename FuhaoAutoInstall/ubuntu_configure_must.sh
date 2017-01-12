#!/bin/bash

#Language support
sudo gnome-language-selector

#Ubuntu-tweak
sudo add-apt-repository ppa:tualatrix/ppa
sudo apt-get update
sudo apt-get install ubuntu-tweak


#NGS library
#Zlib
sudo apt-get install ruby
sudo apt-get install zlib1g zlib1g-dev
#boost
sudo apt-get install libboost-dev libbz2-dev sqlite3 libsqlite3-dev libncurses5-dev libpng-dev libgd-dev libexpat1-dev db6.0
python-dev build-essential gfortran
#svn
sudo apt-get install subversion
#video codes
sudo apt-get install gstreamer0.10-plugins-ugly gxine libdvdread4 totem-mozilla icedax tagtool easytag id3tool lame nautilus-script-audio-convert libmad0 mpg321 gstreamer1.0-libav
#flash plugins
sudo apt-get install flashplugin-installer
#compress and uncompress
sudo apt-get install unace unrar zip unzip p7zip-full p7zip-rar sharutils rar uudeview mpack arj cabextract file-roller

sudo apt-get install dh-autoreconf

#cython python-numpy python-biopython python-biopython-doc

#Chromium
sudo apt-get install chromium-browser

#Gedit-developer-plugins
sudo add-apt-repository ppa:sinzui/ppa
sudo apt-get update
sudo apt-get install gedit-developer-plugins

sudo apt-get install gedit-plugins


#cinnamon 
sudo add-apt-repository ppa:gwendal-lebihan-dev/cinnamon-nightly 
sudo apt-get update 
sudo apt-get install cinnamon



#PDFLatex
#Install the TexLive base 
sudo apt-get install texlive-latex-base
#Also install the recommended and extra fonts to avoid running into the error [1], when trying to use pdflatex on latex files with more fonts.
sudo apt-get install texlive-fonts-recommended
sudo apt-get install texlive-fonts-extra
#Install the extra packages,
sudo apt-get install texlive-latex-extra



#Chinese input
#System Settings\Personal\ Language Support
#	Install/Remove Languages
#	勾选Chinese（simplified），点击“Apply Changes”
#	“Language Support”->“Keyboard input method system” 中选择 ibus
#	点击“Apply System-Wide”
#	System Settings->Personal下的“Text Entry（文本输入）”。
#	点击左下角的“+”，在打开的“Choose an input source”中找到“Chinese(Pinyin)”，
#	点击“Add”，添加。添加后，可以点击“^”“v”调整输入法的默认位置。
sudo apt-get install ibus-googlepinyin
#ibus-setup
#input method, check ”customize active input methods”, select an input method -> show all input methods -> Chinese -> GooglePinyin



#添加中文字符编码的方法
sudo locale-gen zh_CN.GB18030
sudo locale-gen zh_CN.GBK
sudo locale-gen zh_CN.GB2312
dpkg-reconfigure locales
locale-gen

#sudo vim /var/lib/locales/supported.d/zh 加入以下配置参数
#zh_CN.GB18030 GB18030
#zh_CN.GBK GBK
#zh_CN.GB2312 GB2312
#zh_HK.BIG5 BIG5
#zh_TW.BIG5 BIG5
#然后执行 
sudo locale-gen
#提示以下信息，成功了（可能比较慢，耐心等待）
#zh_CN.GB18030... done
#zh_CN.GBK... done



#gedit
gsettings set org.gnome.gedit.preferences.encodings auto-detected "['GB18030', 'GB2312', 'GBK', 'UTF-8', 'BIG5', 'CURRENT', 'UTF-16']"
gsettings set org.gnome.gedit.preferences.encodings shown-in-menu "['GB18030', 'GB2312', 'GBK', 'UTF-8', 'BIG5', 'CURRENT', 'UTF-16']"


#git

sudo apt-get install git





