# Deploy: https://mzdesk.com/

Server 143.198.233.45 (Ubuntu 24.10, nginx in front of MCID / my-trip / mzDesk).

| What | Where |
|---|---|
| Code + built frontend + MetFrag jar | `/opt/mzdesk` (owner root, group mzdesk, read-only for the app) |
| Python venv (core deps, no deep-learning scores, no Docker / ESI MS/MS) | `/opt/mzdesk/.venv` |
| Models, Known-methods indexes, caches, `stats.sqlite` | `/var/lib/mzdesk` (owner mzdesk) |
| Feedback and crash reports (`reports.sqlite`) | `/var/lib/mzdesk/feedback` (owner mzdesk) |
| Secrets and mail settings | `/etc/mzdesk/mzdesk.env` (root, 0600; `mzdesk.env.example`) |
| Spectral libraries (MassBank + GNPS index), CCS data | `/var/lib/mzdesk/data/libraries`, `/var/lib/mzdesk/data/ccs` |
| Engines: Mono (conda-forge, has the Roslyn `csc`) + MS-FINDER, FRIGID (own Python 3.10 env, CPU PyTorch, 2.9 GB of checkpoints) | `/var/lib/mzdesk/tools/` |
| CFM-ID | Docker image `wishartlab/cfmid:msml4.4.10_msrb1.1.15`; the service has the `docker` group (`SupplementaryGroups`) |
| Service | `mzdesk.service` → uvicorn on 127.0.0.1:8599, user `mzdesk`, sandboxed, `Restart=always`, one worker |
| nginx | `/etc/nginx/conf.d/mzdesk.conf` (`nginx-mzdesk.conf`), rate limits in `00-mzdesk-limits.conf`; the old `https://bowenyang.org/mzdesk/` redirects to `mzdesk.com` |

Update after a change: `deploy/deploy.sh` (needs the deploy key `~/.ssh/mzdesk_deploy`). It builds the frontend, syncs the code, installs the
Python packages of `requirements-server.txt` that are missing, installs the service unit and the nginx files when they changed (the old ones go to `/root/nginx-backup/`, `nginx -t` guards the reload), restarts and checks.
The service runs a single uvicorn worker (`--workers 1`): workflow tasks, the queue and the heavy-task lock are kept in process memory, so more
workers would make task polls return 404 until the task store is shared.
Logs: `journalctl -u mzdesk -f`. Restart: `systemctl restart mzdesk`.
The unit sets `MZDESK_HOSTED=1` (this is the public website): `/api/spec` then says `hosted`, the sidebar drops the "Internet search" switch and online lookups
run without asking, since what people type reaches this server anyway. The local version leaves it unset. Read at start, so restart after changing it.

## Releases and the install command

People install the local version with one command (the install page of the Docs shows it):

- macOS / Linux: `curl -fsSL https://mzdesk.com/install.sh | sh`
- Windows (PowerShell): `irm https://mzdesk.com/install.ps1 | iex`

The app serves four things when `MZDESK_HOSTED=1` (`mzdesk/api/release.py`): `/install.sh` and `/install.ps1` (the files in `mzdesk/install/`),
`/version.json`, and `/download/mzdesk-X.Y.Z-py3-none-any.whl` (the wheel, plus its `.sha256` file). The last two are read from
`/var/lib/mzdesk/releases` (`MZDESK_RELEASES_DIR`) at each request, so a release is only a set of files; no restart.
The scripts check the wheel against the `sha256` in `version.json`; `mzdesk update` on a user's machine runs the same script again.

`deploy/build-release.sh` builds `dist/release/` from this checkout (the version is `mzdesk/__init__.py`): the wheel, its checksum and `version.json`.
It needs the built frontend (`cd web && npm run build`) and the MetFrag jar, and it looks inside the wheel before writing anything: the pages, docs,
figures, install scripts, jar and data must be there; benchmark data, tests, secrets and databases must not; large files are listed.
`deploy/deploy.sh` runs it (with the frontend it just built), uploads `dist/release/` to `/root/mzdesk-deploy/release/`, and `remote-install.sh`
publishes it **after** the site answers: the wheel and checksum first, then `version.json` in one `mv`, so the old release stays whole until the new one is complete.
The folder `/var/lib/mzdesk/releases` belongs to root (group `mzdesk`, read only for the app), and the newest five wheels are kept, so
`curl -fsSL https://mzdesk.com/install.sh | MZDESK_VERSION=0.5.3 sh` still installs an older one.

A release is a version bump: change `__version__`, deploy. Run `deploy/build-release.sh` alone to test the wheel (nothing is uploaded). With
`MZDESK_RELEASE_ALLOW_STUB_STATIC=1` it builds without a frontend, for tests of the installer only; `deploy.sh` refuses that setting.

nginx: `/download/` has its own limit zone (`mzdesk_download`, 6 per minute per address, burst 6), no response buffering and 15 minutes between two
writes to a slow client; the install scripts and `version.json` fall under the normal page limit. The app sends the wheel with `Content-Encoding: identity`
so that its gzip middleware leaves the already compressed file alone. After the first deploy, check from another machine:
`curl -fsSL https://mzdesk.com/install.sh | head -3` and `curl -s https://mzdesk.com/version.json`.

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

## The engines (all pages run here, none is "local only")

`deploy/install-engines.sh` installs what the local version brings with it (about 7 GB downloaded, 10 GB on disk; FRIGID's checkpoint from Zenodo takes the longest);
every step skips what is in place. The Python packages of the main venv (pyopenms, PyTorch CPU, ms2deepscore, gensim, ...) come with `deploy.sh`
(`requirements-server.txt`). The service unit points mzDesk at all of it (`MZDESK_*_DIR`, `MZDESK_MONO`).

- **FRIGID** runs on the CPU: a Quick run takes about two minutes and peaks at **3.1 GB resident** (measured), on a 3.9 GB server. The unit therefore allows
  the service 3.4 GB (`MemoryMax`) and sets `OOMScoreAdjust=300`, so if memory ever runs out the kernel stops mzDesk first and not MCID or the other sites.
  During a run the rest of the box is squeezed (swap reached 2.1 GB); the other sites kept answering in 0.2 to 0.5 s. One heavy job (FRIGID or MS-FINDER) runs at
  a time. MS/MS Annotation's "Run FRIGID" (Advanced) starts unticked on the server (`MZDESK_ANNOTATION_FRIGID=off` in the unit) and ticked in the local version;
  the De Novo Structure page always runs it when asked. More RAM on the droplet is the real fix. nginx lets the De Novo page wait 1300 s (the engine stops a Full run at 20 minutes).
- **CFM-ID** runs in Docker. The service may start containers, which is root-equivalent if the app is ever taken over, so the app never builds a Docker command from user
  text: the image and arguments are fixed in `ms2/cfmid_engine.py`, the container gets no network, no capabilities, no privilege gain, a process cap, 1 GB of memory
  and runs as the service's own user. The run's folder is `MZDESK_CFMID_WORKDIR` because the unit's `/tmp` is private and the Docker daemon cannot see it.
- Ubuntu's Mono package has no `csc`, which MS-FINDER's console front end needs; Mono comes from conda-forge into `/var/lib/mzdesk/tools/mono`.
- The deep-learning scores download their models on first use (MS2DeepScore 416 MB, Spec2Vec about 280 MB, from Zenodo, a few minutes): the first request that
  asks for them can time out in nginx while the download finishes in the background; the next one is fast.

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
should have its own: `deploy/set-admin-password.sh` asks for a password (hidden, twice), puts only its salted hash into `/etc/mzdesk/mzdesk.env`, restarts the service
and checks `/stats` and `/reports` with the new and a wrong password. (`python -m mzdesk.api.admin` only prints the two `MZDESK_ADMIN_*` lines.)

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
