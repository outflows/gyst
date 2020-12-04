#!/bin/bash

for i in $(seq -f "%03g" 0 600)
do
    ./gyst "dump$i"
done
