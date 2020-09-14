#!/bin/bash

for i in $(seq -f "%03g" 0 700)
do
    ./gyst "dump$i"
done
