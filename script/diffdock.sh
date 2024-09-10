# docker run -d --gpus all -p 7860:7860 --entrypoint /bin/bash rbgcsail/diffdock -c "micromamba run -n diffdock python app/main.py"

# # 启动diffdock在线
# docker run -d  -p 7860:7860 --entrypoint /bin/bash dockerhub.icu/rbgcsail/diffdock:latest -c "micromamba run -n diffdock python app/main.py"

# # 始终重启
# docker run -d --restart unless-stopped -p 7860:7860 --entrypoint /bin/bash dockerhub.icu/rbgcsail/diffdock:latest -c "micromamba run -n diffdock python app/main.py"

# 固定模型的版本，免得临时下载模型失败。 等等成功
docker run -d --restart unless-stopped -p 7860:7860 --entrypoint /bin/bash august777/diffdocker:full -c "micromamba run -n diffdock python app/main.py"
