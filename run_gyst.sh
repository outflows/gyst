#!/bin/bash

for i in $(seq -f "%03g" 0 1000)
do
    ./gyst "dump$i"
done
