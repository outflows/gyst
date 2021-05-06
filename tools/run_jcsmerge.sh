#!/bin/bash

for i in $(seq -f "%03g" 225 940)
do
    ./jcsmerge "dump$i"
done
