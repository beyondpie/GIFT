#!/bin/bash

# Test package gift.
# Songpeng Zu
# 2016-03-04

# set executable gift
gift=../bin/gift
# test help func
${gift} -h
${gift} -help
# test gift version
${gift} -version
${gift} -v
# test gift train
${gift} -config test_init-train.gift -d2p test_drug2protein \
        -d2s test_drug2sub -p2d test_protein2domain
# test gift predict, drug
cp chemFP2proFP* ./test_s2d
${gift} -config test_init-predict.gift -s2d test_s2d \
        -d2s test_drug2sub -druglist test_drugNameList
# test gift predict, protein
${gift} -config test_init-predict.gift -s2d test_s2d \
        -p2d test_protein2domain -proteinlist test_proteinNameList
