#!/bin/bash

while (( 1 < 2 ))
do

    echo "Starting scan"
    timestamp=$(date +%s)
    sudo timeout 3600s tcpdump -i enp0s31f6 -w $timestamp.pcapng udp port 6006
    echo "Finished scan"
    sleep 5
done

