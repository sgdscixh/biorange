# 安装构建工具
conda install anaconda-client conda-build 

# 安装conda 最先进得打包工具
pip install grayskull 

# 直接转换poetry得包（很轻松几下就结束了，但是需要处理兼容问题）
grayskull pypi ./dist/biorange-1.4.1.tar.gz 

conda-build  biorange   

conda build purge 
anaconda upload  biorange-1.3.0-py_0.tar.bz2 