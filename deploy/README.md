# Deploy: https://bowenyang.org/mzdesk/

Server 143.198.233.45 (Ubuntu 24.10, nginx in front of MCID / my-trip / mzDesk).

| What | Where |
|---|---|
| Code + built frontend + MetFrag jar | `/opt/mzdesk` (owner root, group mzdesk, read-only for the app) |
| Python venv (core deps, no deep-learning scores, no Docker / ESI MS/MS) | `/opt/mzdesk/.venv` |
| Models, Known-methods indexes, caches, `stats.sqlite` | `/var/lib/mzdesk` (owner mzdesk) |
| Service | `mzdesk.service` → uvicorn on 127.0.0.1:8599, user `mzdesk`, sandboxed, `Restart=always` |
| nginx | `location /mzdesk/` in `/etc/nginx/conf.d/bowenyang.conf` (`nginx-location.conf`), rate limits in `00-mzdesk-limits.conf` |

Update after a change: `deploy/deploy.sh` (needs the deploy key `~/.ssh/mzdesk_deploy`).
Logs: `journalctl -u mzdesk -f`. Restart: `systemctl restart mzdesk`. nginx backups: `/root/nginx-backup/`.

Security: app listens on localhost only; per-IP rate limits (pages 20/s, module runs 30/min, stats password 6/min);
8 MB request cap; nosniff / frame / referrer / permissions headers; `/docs` and `/openapi.json` hidden; no IPs stored in stats.
The page lock is still frontend-only for pages; the stats numbers are protected on the server with the same password.
