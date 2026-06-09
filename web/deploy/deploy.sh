#!/bin/bash
# ============================================================
# UVTOSIZE Deployment Script
# Run on the liuzunyu.cn server as root or with sudo.
#
# Usage:
#   chmod +x deploy.sh
#   sudo ./deploy.sh
# ============================================================
set -e

echo "=============================="
echo " UVTOSIZE Deployment"
echo "=============================="

# ---- Configuration ----
APP_DIR="/var/www/uvtosize"
VENV_DIR="$APP_DIR/venv"
STATIC_DIR="/var/www/liuzunyu/web/static"
NGINX_CONF="/etc/nginx/sites-available"
LOG_DIR="/var/log/uvtosize"
PORT=8765

# ---- 1. Create directories ----
echo "[1/7] Creating directories..."
mkdir -p "$APP_DIR"
mkdir -p "$APP_DIR/scripts"
mkdir -p "$STATIC_DIR"
mkdir -p "$LOG_DIR"
chown -R www-data:www-data "$LOG_DIR"

# ---- 2. Copy application files ----
echo "[2/7] Copying application files..."
# Copy server and WSGI entry point
cp server.py "$APP_DIR/"
cp wsgi.py "$APP_DIR/"
cp requirements_web.txt "$APP_DIR/"

# Copy the analysis script
cp ../.claude/skills/UVTOSIZE/scripts/uv_analysis.py "$APP_DIR/scripts/"

# Copy the frontend HTML
cp static/uvtosize.html "$STATIC_DIR/"

# ---- 3. Set up Python virtual environment ----
echo "[3/7] Setting up Python virtual environment..."
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip
pip install -r "$APP_DIR/requirements_web.txt"
pip install gunicorn

# ---- 4. Set permissions ----
echo "[4/7] Setting permissions..."
chown -R www-data:www-data "$APP_DIR"
chmod -R 755 "$APP_DIR"

# ---- 5. Configure systemd service ----
echo "[5/7] Installing systemd service..."
cp deploy/uvtosize.service /etc/systemd/system/
systemctl daemon-reload
systemctl enable uvtosize.service

# ---- 6. Configure Nginx ----
echo "[6/7] Configuring Nginx..."
echo ""
echo "  >> Add the following to your Nginx server block for liuzunyu.cn:"
echo ""
echo "  location = /uvtosize {"
echo "      alias $STATIC_DIR/uvtosize.html;"
echo "  }"
echo ""
echo "  location /api/uvtosize/ {"
echo "      proxy_pass http://127.0.0.1:$PORT/api/uvtosize/;"
echo "      proxy_http_version 1.1;"
echo "      proxy_set_header Host \$host;"
echo "      proxy_set_header X-Real-IP \$remote_addr;"
echo "      proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;"
echo "      proxy_set_header X-Forwarded-Proto \$scheme;"
echo "      proxy_read_timeout 120s;"
echo "      proxy_send_timeout 120s;"
echo "      client_max_body_size 50m;"
echo "  }"
echo ""
echo "  Full config at: deploy/nginx-uvtosize.conf"
echo "  Then run: sudo nginx -t && sudo systemctl reload nginx"

# ---- 7. Start service ----
echo "[7/7] Starting UVTOSIZE service..."
systemctl restart uvtosize.service
sleep 2
systemctl status uvtosize.service --no-pager

echo ""
echo "=============================="
echo " Deployment Complete!"
echo "=============================="
echo ""
echo "  Tool page:  https://liuzunyu.cn/uvtosize"
echo "  API health: https://liuzunyu.cn/api/uvtosize/health"
echo "  Service:    sudo systemctl [status|restart|stop] uvtosize"
echo "  Logs:       sudo journalctl -u uvtosize -f"
echo "              tail -f $LOG_DIR/access.log"
echo ""
