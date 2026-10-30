# Deploy: https://bowenyang.org/mzdesk/

Server 143.198.233.45 (Ubuntu 24.10, nginx in front of MCID / my-trip / mzDesk).

| What | Where |
|---|---|
| Code + built frontend + MetFrag jar | `/opt/mzdesk` (owner root, group mzdesk, read-only for the app) |
| Python venv (core deps, no deep-learning scores, no Docker / ESI MS/MS) | `/opt/mzdesk/.venv` |
| Models, Known-methods indexes, caches, `stats.sqlite` | `/var/lib/mzdesk` (owner mzdesk) |
| Feedback awaiting email delivery (`.eml`) | `/var/lib/mzdesk/feedback` (owner mzdesk) |
| Service | `mzdesk.service` → uvicorn on 127.0.0.1:8599, user `mzdesk`, sandboxed, `Restart=always` |
| nginx | `location /mzdesk/` in `/etc/nginx/conf.d/bowenyang.conf` (`nginx-location.conf`), rate limits in `00-mzdesk-limits.conf` |

Update after a change: `deploy/deploy.sh` (needs the deploy key `~/.ssh/mzdesk_deploy`).
Logs: `journalctl -u mzdesk -f`. Restart: `systemctl restart mzdesk`. nginx backups: `/root/nginx-backup/`.

GNPS2: install the updated `requirements.txt` in the server venv before restarting (MassQL mzML/mzXML needs
`pyopenms`). Apply the updated `nginx-location.conf` and reload nginx: only MassQL's endpoint allows the 20 MB file
(28 MiB JSON request cap); fastMASST has a longer read timeout. The 2026-10-12 Pan-ReDU snapshot ships with the code;
refreshed caches go under `/var/lib/mzdesk/cache/gnps/`. Details and local CLI: `admin/GNPS2.md`.

## Feedback email

The sidebar's bug button posts to `/api/feedback`; the server emails the message and up to three screenshots to
`MZDESK_FEEDBACK_TO` (default `feedback@example.org`). SMTP credentials live in `/etc/mzdesk/mail.env`, which the unit
loads with `EnvironmentFile=-/etc/mzdesk/mail.env`:

```bash
install -d -m 750 /etc/mzdesk
install -m 600 deploy/mail.env.example /etc/mzdesk/mail.env
$EDITOR /etc/mzdesk/mail.env                  # fill in host / user / password / from
systemctl daemon-reload && systemctl restart mzdesk
```

Use any authenticated SMTP relay. A transactional provider (Resend, Postmark, SendGrid, Mailgun) or a Gmail app
password is usually easier than Georgetown's Microsoft 365, which may require OAuth and disable basic SMTP AUTH.
DigitalOcean permits outbound 587. `MZDESK_SMTP_FROM` must be an address the relay is allowed to send as.

If SMTP is missing or fails, nothing is lost: the complete message is written as an `.eml` file under
`/var/lib/mzdesk/feedback/`, the form still says "saved", and the failure is logged (`journalctl -u mzdesk | grep feedback`).

The app also rate-limits feedback (30 s between messages, 8/hour per IP per worker) and has a hidden honeypot field.
The nginx rule (`mzdesk_feedback`, 6/min) is an additional shared limit; reload nginx after changing
`nginx-location.conf` / `00-mzdesk-limits.conf`:

```bash
nginx -t && systemctl reload nginx
```

Security: app listens on localhost only; per-IP rate limits (pages 20/s, module runs 30/min, stats password 6/min,
feedback 6/min); 8 MB request cap; nosniff / frame / referrer / permissions headers; `/docs` and `/openapi.json`
hidden; no IPs stored in stats; feedback images are validated by their bytes and only ever attached to the email.
The page lock is still frontend-only for pages; the stats numbers are protected on the server with the same password.
