# Deploy: https://mzdesk.com/

Server 143.198.233.45 (Ubuntu 24.10, nginx in front of MCID / my-trip / mzDesk).

| What | Where |
|---|---|
| Code + built frontend + MetFrag jar | `/opt/mzdesk` (owner root, group mzdesk, read-only for the app) |
| Python venv (core deps, no deep-learning scores, no Docker / ESI MS/MS) | `/opt/mzdesk/.venv` |
| Models, Known-methods indexes, caches, `stats.sqlite` | `/var/lib/mzdesk` (owner mzdesk) |
| Feedback and crash reports (`reports.sqlite`) | `/var/lib/mzdesk/feedback` (owner mzdesk) |
| Secrets and mail settings | `/etc/mzdesk/mzdesk.env` (root, 0600; `mzdesk.env.example`) |
| Service | `mzdesk.service` → uvicorn on 127.0.0.1:8599, user `mzdesk`, sandboxed, `Restart=always`, one worker |
| nginx | `/etc/nginx/conf.d/mzdesk.conf` (`nginx-mzdesk.conf`), rate limits in `00-mzdesk-limits.conf`; the old `https://bowenyang.org/mzdesk/` redirects to `mzdesk.com` |

Update after a change: `deploy/deploy.sh` (needs the deploy key `~/.ssh/mzdesk_deploy`). It builds the frontend, syncs the code, installs the
Python packages of `requirements-server.txt` that are missing, installs the service unit and the nginx files when they changed (the old ones go to `/root/nginx-backup/`, `nginx -t` guards the reload), restarts and checks.
The service runs a single uvicorn worker (`--workers 1`): workflow tasks, the queue and the heavy-task lock are kept in process memory, so more
workers would make task polls return 404 until the task store is shared.
Logs: `journalctl -u mzdesk -f`. Restart: `systemctl restart mzdesk`.

## The address (done once)

mzdesk.com is registered at Porkbun. At Porkbun, under DNS for mzdesk.com, the records are:

| Type | Host | Answer |
|---|---|---|
| A | (empty) | 143.198.233.45 |
| A | www | 143.198.233.45 |

Remove Porkbun's default parking records (an ALIAS / CNAME on the bare name and on `www`), and leave the MX records alone (they forward mail).
Then run `deploy/cutover.sh` from this checkout. It checks the DNS, gets the Let's Encrypt certificate (webroot, renewed by the existing
`certbot.timer`), deploys, installs the https site and turns `bowenyang.org/mzdesk/...` into a 301 to `mzdesk.com/...`.

The frontend is built for the root (`/`). `MZDESK_BASE=/sub/` still builds it for a sub-path.

## Feedback and crash reports

Both go through `mzdesk/api/reports.py`:

- **Feedback**: the sidebar's bug button posts to `/api/feedback` (message, up to three screenshots).
- **Crashes**: a run that raises on the server (module, workflow step, fragment view, any `log.exception` in mzdesk) and a page that fails to draw in
  the browser (`/api/crash`) are recorded with the traceback. The user sees a reference (`E-1A2B3C`) and nothing of the server's folders. The same crash
  within a day counts up the first report and is forwarded at most once an hour. Of the inputs of a run only the field names and sizes are kept.

Every report is **stored first** in `/var/lib/mzdesk/feedback/reports.sqlite` (the newest 500 are kept) and read at **https://mzdesk.com/reports**
(admin password; delete a report there when it is dealt with). This needs no setup and cannot fail because a mail service is down.

Forwarding is optional and independent per channel; set it in `/etc/mzdesk/mzdesk.env` and `systemctl restart mzdesk`:

- **Email through Resend** (`MZDESK_RESEND_API_KEY`): this droplet cannot send by SMTP (DigitalOcean blocks outbound ports 25, 465 and 587), so email goes over
  HTTPS. Sign up at resend.com with the address that should receive the reports, create an API key, put it in the file. Without a verified domain Resend
  delivers only to that address, from `onboarding@resend.dev`, which is what a report to yourself needs. To send from `reports@mzdesk.com`, verify mzdesk.com in
  Resend (it shows a few DNS records to add at Porkbun) and set `MZDESK_MAIL_FROM`.
- **Push through ntfy** (`MZDESK_NTFY_URL`): no account; choose a long random topic, subscribe to it in the ntfy app. The push carries the reference, the kind and
  the page, never the message text.
- **SMTP**: `MZDESK_SMTP_*`, for a server that can reach one.

## The admin password

`/stats` and `/reports` check the password on the server, with a limit of 8 wrong guesses per address in 15 minutes (nginx limits the paths too).
By default that is the page lock's password, whose hash ships to every browser inside the page, so a short password can be guessed offline. A server
should have its own: `python -m mzdesk.api.admin` asks for a password and prints the two `MZDESK_ADMIN_*` lines for `/etc/mzdesk/mzdesk.env`.

## Security

- The app listens on localhost only; nginx is the only way in. TLS 1.2 / 1.3 (certbot's options), HSTS on mzdesk.com, `server_tokens off` there.
- Headers on every response (`mzdesk/api/security.py`): `Content-Security-Policy` (own origin only, scripts by hash, nothing may frame the app),
  `X-Content-Type-Options`, `X-Frame-Options`, `Referrer-Policy`, `Permissions-Policy`, `Cross-Origin-Opener-Policy`; API answers are `no-store`.
- Request bodies are capped in the app (8 MB, MassQL 27 MB) as well as in nginx; per-address rate limits (pages 20/s, module runs 30/min, admin password 6/min,
  feedback 6/min, crash reports 30/min); the app limits feedback (30 s between messages, 8 an hour) and crash reports (30 an hour) too.
- Errors shown to people carry no folders of the server; the traceback stays in the log and in `/reports`.
- FastAPI's `/docs`, `/redoc` and `/openapi.json` are off (`MZDESK_API_DOCS=1` turns them on for local development).
- No IPs in stats; screenshots are checked by their bytes and stored only for the maintainer.
- The service unit is sandboxed: no new privileges, read-only system, one writable folder, no devices, no kernel access, no other users' processes.
- Dependencies were audited (pip-audit on the server's packages: nothing known). `npm audit` lists react-router 6 (open redirect through
  user-controlled `<Link>` targets and server-side hydration, neither used) and Ketcher's build-time glob; their fixes are major upgrades.

GNPS2: install the updated `requirements.txt` in the server venv before restarting (MassQL mzML/mzXML needs `pyopenms`). Only MassQL's endpoint allows
the 20 MB file; fastMASST has a longer read timeout. The 2026-10-12 Pan-ReDU snapshot ships with the code; refreshed caches go under
`/var/lib/mzdesk/cache/gnps/`. Details and local CLI: `admin/GNPS2.md`.
