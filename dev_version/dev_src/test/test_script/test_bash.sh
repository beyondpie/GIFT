#!/bin/bash

# Test package gift.
# Songpeng Zu
# 2016-03-04

# set executable gift
gift=/Users/wangchao/home/songpeng/git-recipes/GIFT/dev_version/dev_src/build/bin/gift
# test help func
${gift} -h
${gift} -help
# test gift version
${gift} -version
${gift} -v

# test gift train
${gift} --configure test_init-train.gift
