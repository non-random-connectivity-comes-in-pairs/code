#!/bin/bash

docker run -it -p 127.0.0.1:8017:8000 \
       --user="$(id -u):$(id -g)" \
       -v $(pwd):/home/lab \
       felix11h/scipy3_env \
       /bin/bash -c \
       'source startup_messg.sh;
        bash'
