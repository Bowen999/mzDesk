"""The local version of mzDesk (the `mzdesk` command on a user's own computer): where its engines live, which first-use downloads
("components") exist, what each page needs, and how a download runs. Nothing here is used on mzdesk.com (MZDESK_HOSTED=1).

  paths.py       folders (~/.mzdesk/...), Java and git discovery; standard library only, so the engines can import it
  archive.py     files of one pinned GitHub commit without git (CCS and RepoRT setups)
  components.py  the registry: id, size, installed(), install()
  targets.py     which page / workflow needs which component
  prereq.py      what an install provisions first: micromamba, Mono, git, a throw-away uv environment, the Docker check
  runner.py      child processes: stream their output, kill the tree on cancel
  jobs.py        one install job at a time, progress and log tail

The package is imported by the engines (paths) and by the API (everything else); keep this file free of imports.
"""
