docker run \
   -d -p 6666:9000 \
   -p 9011:9001 \
   --name minio \
   -v /home/liuyan/database/:/data \
   -e "MINIO_ROOT_USER=liuyan" \
   -e "MINIO_ROOT_PASSWORD=liuyan666" \
   quay.io/minio/minio server /data --console-address ":9001"