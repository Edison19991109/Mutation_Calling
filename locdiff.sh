#!/usr/bin/env bash
set -ue
[ -f "${1}" ] && [ -f "${2}" ]
echo "Item in ${1} but not ${2}:"
cat "${1}".noend "${2}".noend "${2}".noend |sort|uniq -u
echo "Item in ${2} but not ${1}:"
cat "${2}".noend "${1}".noend "${1}".noend |sort|uniq -u
