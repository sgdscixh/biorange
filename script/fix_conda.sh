mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh 

# conda 坏了使用这个命令修复，目前坏得原因是用conda在base环境升级conda，包来源不同一个是anaconda 一个是conda-forge自然崩溃了 后面他们会修复
