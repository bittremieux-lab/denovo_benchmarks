#!/bin/bash
SERV_PID=$(ps -ef|grep novorai_server|grep python|awk '{print $2}')
[ -n "${SERV_PID}" ] && kill -TERM ${SERV_PID}
