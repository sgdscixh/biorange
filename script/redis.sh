docker run --name redis  -d \
  --restart unless-stopped \
  -p 6379:6379 \
  -v redis-data:/data \
  redis:latest
# 自建的docker 镜像
# https://docker.fanxingfu3344.workers.dev/