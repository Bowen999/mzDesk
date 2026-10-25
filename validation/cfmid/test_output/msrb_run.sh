#!/bin/sh
cd /data/timing || exit 1
java -jar /cfmid/msrb-fragmenter.jar -a "$1" -ismi "$2" -o "$3" 2>&1
