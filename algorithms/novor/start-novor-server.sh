#!/bin/bash
LOGFILE="novor-gnn-server.log"
CHKRUN="$(ps -ef|grep novorai_server|grep python)"
[ -n "${CHKRUN}" ] && exit 0
nohup /home/novor/bin/novor-gnn-server > $LOGFILE 2>&1 &
while true; do
    sleep 1
    [ -f "$LOGFILE" ] || break
    [ -n "$(grep -o 'Server starts.' $LOGFILE)" ] && break
done
echo "Novor server started."