# UVTOSIZE 生产部署指南

将量子点 UV-Vis 分析工具部署到 liuzunyu.cn 服务器。

## 架构

```
用户浏览器                    Nginx                       Flask (Gunicorn)
    │                          │                             │
    ├─ /uvtosize ──────────────┤─ static HTML ───────────────┤
    ├─ /api/uvtosize/analyze ──┤─ reverse proxy ─────────────┤─ analyze_to_dict()
    ├─ /api/uvtosize/download ─┤─ reverse proxy ─────────────┤─ send_file()
    └─ /tools (existing) ──────┤─ static pages ──────────────┤
```

## 目录结构 (服务器上)

```
/var/www/
├── liuzunyu/                     # 现有站点
│   └── web/static/
│       └── uvtosize.html         # ← 从 web/static/ 上传
└── uvtosize/                     # 新建
    ├── server.py                 # Flask 应用
    ├── wsgi.py                   # Gunicorn 入口
    ├── requirements_web.txt      # Python 依赖
    ├── scripts/
    │   └── uv_analysis.py        # 分析引擎
    ├── venv/                     # Python 虚拟环境
    └── deploy/
        ├── nginx-uvtosize.conf
        ├── uvtosize.service
        └── deploy.sh

/var/log/uvtosize/                # 日志
    ├── access.log
    └── error.log
```

## 快速部署 (推荐)

### 1. 上传文件到服务器

在本地执行 (从 `web/` 目录):

```bash
# 上传应用文件
scp -r server.py wsgi.py requirements_web.txt user@liuzunyu.cn:/var/www/uvtosize/

# 上传分析脚本
scp ../.claude/skills/UVTOSIZE/scripts/uv_analysis.py \
    user@liuzunyu.cn:/var/www/uvtosize/scripts/

# 上传前端页面
scp static/uvtosize.html user@liuzunyu.cn:/var/www/liuzunyu/web/static/

# 上传部署配置
scp -r deploy/ user@liuzunyu.cn:/var/www/uvtosize/deploy/
```

### 2. SSH 到服务器并运行部署脚本

```bash
ssh user@liuzunyu.cn
sudo bash /var/www/uvtosize/deploy/deploy.sh
```

### 3. 配置 Nginx

将以下内容添加到站点 Nginx 配置中 (通常在 `/etc/nginx/sites-available/liuzunyu`):

```nginx
# === UVTOSIZE 工具页面 ===
location = /uvtosize {
    alias /var/www/liuzunyu/web/static/uvtosize.html;
}

# === UVTOSIZE API ===
location /api/uvtosize/ {
    proxy_pass http://127.0.0.1:8765/api/uvtosize/;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_read_timeout 120s;
    proxy_send_timeout 120s;
    client_max_body_size 50m;
}
```

然后:

```bash
sudo nginx -t
sudo systemctl reload nginx
```

### 4. 验证

```bash
# 检查服务状态
sudo systemctl status uvtosize

# 测试 API
curl https://liuzunyu.cn/api/uvtosize/health

# 在浏览器打开
# https://liuzunyu.cn/uvtosize
```

## 手动部署步骤

### 前置条件

- Ubuntu 20.04+ / Debian 11+ (服务器)
- Python 3.10+
- Nginx
- systemd
- `www-data` 用户

### 详细步骤

#### 1. 创建目录和用户

```bash
sudo mkdir -p /var/www/uvtosize/scripts
sudo mkdir -p /var/www/uvtosize/deploy
sudo mkdir -p /var/log/uvtosize
sudo chown -R www-data:www-data /var/www/uvtosize /var/log/uvtosize
```

#### 2. 复制文件

```bash
# 在项目目录本地执行
rsync -avz --progress \
  server.py wsgi.py requirements_web.txt \
  user@liuzunyu.cn:/var/www/uvtosize/

rsync -avz --progress \
  ../.claude/skills/UVTOSIZE/scripts/uv_analysis.py \
  user@liuzunyu.cn:/var/www/uvtosize/scripts/

rsync -avz --progress \
  static/uvtosize.html \
  user@liuzunyu.cn:/var/www/liuzunyu/web/static/
```

#### 3. 安装 Python 依赖

```bash
ssh user@liuzunyu.cn
cd /var/www/uvtosize
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements_web.txt
pip install gunicorn
```

#### 4. 测试应用

```bash
cd /var/www/uvtosize
source venv/bin/activate
gunicorn wsgi:app --bind 127.0.0.1:8765 --workers 2
# Ctrl+C 停止
```

#### 5. 安装 systemd 服务

```bash
sudo cp /var/www/uvtosize/deploy/uvtosize.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable uvtosize
sudo systemctl start uvtosize
sudo systemctl status uvtosize
```

#### 6. 配置 Nginx (见上方)

#### 7. 在工具中心添加入口

在 `/var/www/liuzunyu/web/static/tools.html` 或对应页面，在工具卡片区域添加:

```html
<a href="/uvtosize" class="tool-card">
  <div class="tool-icon">🔬</div>
  <div class="tool-name">UVTOSIZE · 量子点分析</div>
  <div class="tool-desc">从UV-Vis吸收光谱自动计算量子点尺寸与粒径分布</div>
  <div class="tool-price">免费</div>
</a>
```

并在工具切换栏添加:

```html
<a href="/uvtosize">UVTOSIZE</a>
```

## 运维命令

```bash
# 服务管理
sudo systemctl start uvtosize       # 启动
sudo systemctl stop uvtosize        # 停止
sudo systemctl restart uvtosize     # 重启
sudo systemctl status uvtosize      # 状态

# 日志
sudo journalctl -u uvtosize -f      # 实时日志
sudo journalctl -u uvtosize -n 50   # 最近 50 条
tail -f /var/log/uvtosize/access.log
tail -f /var/log/uvtosize/error.log

# 更新代码后重启
sudo systemctl restart uvtosize
sudo systemctl reload nginx
```

## 环境变量

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `PORT` | 8765 | Flask 监听端口 |
| `FLASK_DEBUG` | 0 | 调试模式 (生产环境务必设为 0) |
| `UVTOSIZE_AUTH_REQUIRED` | 0 | 是否需要登录 (1=需要, 0=不需要) |
| `CORS_ORIGINS` | * | 允许的跨域来源 (逗号分隔) |

生产环境建议:

```bash
# 在 /etc/systemd/system/uvtosize.service 的 [Service] 段设置:
Environment="PORT=8765"
Environment="UVTOSIZE_AUTH_REQUIRED=0"
Environment="CORS_ORIGINS=https://liuzunyu.cn"
```

## 故障排查

**服务启动失败**: `sudo journalctl -u uvtosize -n 30`

**502 Bad Gateway**: Flask 未运行 → `sudo systemctl status uvtosize`

**文件上传失败**: 检查 Nginx `client_max_body_size` 是否足够大

**分析报错**: 查看 `/var/log/uvtosize/error.log`

**端口占用**: `sudo lsof -i :8765` → 终止冲突进程
