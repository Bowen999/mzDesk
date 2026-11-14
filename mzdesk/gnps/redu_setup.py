"""Pan-ReDU metadata, one-time download: python -m mzdesk.gnps.redu_setup [--check] [--force]

Dataset Finder (and fastMASST with "Join Pan-ReDU") read the Pan-ReDU per-file metadata of GNPS2. mzDesk does not ship it: GNPS2 states no
licence for redistributing it. This streams the full table from https://redu.gnps2.org/dump (about 1.7 GB, updated daily; a dropped connection is resumed), keeps the 16
columns the pages use and writes one gzip file of about 8 MB (MZDESK_GNPS_REDU_CACHE, else ~/.mzdesk/cache/gnps/redu_slim.tsv.gz).
The full table is never written to disk. --force downloads it again (a newer snapshot); --check only reports.
"""
import argparse
import sys

from mzdesk.gnps import compound_db_check as g


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m mzdesk.gnps.redu_setup", description="Download the Pan-ReDU metadata used by Dataset Finder.")
    ap.add_argument("--check", action="store_true", help="only report whether the metadata is there")
    ap.add_argument("--force", action="store_true", help="download it again even when it is there")
    args = ap.parse_args(argv)
    if args.check:
        if g.redu_installed():
            print(f"Pan-ReDU metadata: {g.REDU_SLIM} (snapshot of {g.redu_cache_date(g.REDU_SLIM)})")
            return 0
        print(f"Pan-ReDU metadata: NOT FOUND ({g.REDU_SLIM})")
        return 1
    if g.redu_installed() and not args.force:
        print(f"Pan-ReDU metadata is already there: {g.REDU_SLIM} (snapshot of {g.redu_cache_date(g.REDU_SLIM)}). --force downloads it again.")
        return 0
    g.redu_download()
    return 0 if g.redu_installed() else 1


if __name__ == "__main__":
    sys.exit(main())
