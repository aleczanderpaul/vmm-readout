#!/bin/bash

echo "Starting scan"
timestamp=$(date +%s)
sudo timeout 1800s tcpdump -i enp0s31f6 -w $timestamp.pcapng udp port 6006
echo "Finished scan"

