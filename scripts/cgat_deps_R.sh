#!/usr/bin/env bash
#
# Script to find R deps in this repository
#
# It simply looks for the 'library(' statement in the source code
#

SCRIPT_FOLDER=$(dirname $0)
REPO_FOLDER=$(dirname ${SCRIPT_FOLDER})

grep -i 'library(' -r ${REPO_FOLDER}/{CGAT,R} \
 | grep -v Binary \
 | sed -e 's/\(.*\)library\(.*\)/\2/' \
 | sed 's/[()"&,.%'\'']//g' \
 | sed 's/\\n$//g' \
 | egrep '^[a-zA-Z]{2,}' \
 | egrep -v 'spp|zinba' \
 | sort -u

