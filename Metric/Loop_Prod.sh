#!/bin/bash

for seas in {1..10}
do
python metrics.py --season ${seas} --type all
done