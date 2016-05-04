#!/bin/bash

# Test package gift.
# Songpeng Zu
# 2016-03-04

pwd
# set executable gift
gift=/Users/wangchao/home/songpeng/git-recipes/GIFT/dev/build/bin

# Compile gift.
gift_dir=/Users/wangchao/home/songpeng/git-recipes/GIFT/dev/
gift_dir_build=${gift_dir}/build
cd ${gift_dir_build}
cmake ..
make

# test help func
# ${gift} -h
# ${gift} --help
# test gift version
# ${gift} --version
# ${gift} -v

# test gift train
cd -
pwd
${gift}/gift --configure test_init-train.gift
