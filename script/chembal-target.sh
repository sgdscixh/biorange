# 原神启动
# docker run -p 8080:8080 a.ussh.net/chembl/mcp

# 本地运行chembal 靶点预测
docker run  --restart unless-stopped -d -p 8080:8080 -p 8081:8081 biotree/mcp:latest
# -d 后台运行 -p 转发端口
# 你说~~我想着先把所有成分都预测了 没问题的
# 测试
curl -X POST -H 'Accept: */*' -H 'Content-Type: application/json' -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}' http://127.0.0.1:8080/predict