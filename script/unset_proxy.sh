# 取消代理
unset https_proxy
unset http_proxy
unset all_proxy  
unset HTTPS_PROXY
unset HTTP_PROXY
unset ALL_PROXY

# 临时在命令行设置代理
export https_proxy=http://127.0.0.1:7890 http_proxy=http://127.0.0.1:7890 all_proxy=socks5://127.0.0.1:7890