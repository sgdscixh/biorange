# 安装构建工具
conda install anaconda-client conda-build 

# 安装conda 最先进得打包工具
pip install grayskull 


# 直接转换poetry得包（很轻松几下就结束了，但是需要处理兼容问题）
poetry build

grayskull pypi ./dist/biorange-1.4.2.tar.gz 
# 这一步反复错没问题的 毕竟两个平台 慢慢解决只有

conda-build  biorange   

conda build purge 
# 本地安装 先测试 本地安装太惨了 不安装依赖
conda install --use-local biorange  --force-reinstall # 通过名字安装
conda install  path/to/xxx.tar.bz2   --force-reinstall# 通过地址安装

# 测试完整再上传

 anaconda upload  ~/miniconda3/conda-bld/noarch/biorange-1.4.2-pypy_0.tar.bz2 --force
