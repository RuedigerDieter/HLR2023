#!/bin/bash
X=$(hostname --short)
Y=$(date --iso-8601=ns)
echo "$X:$Y"
