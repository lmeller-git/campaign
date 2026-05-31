#!/bin/bash


uv run validate.py --keep-work --end 50

cp -r ./metrics /mnt/s3
cp -r ./work /mnt/s3
