# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:41:27 2018

@author: annabelbeichman
"""


def findAllHetLines(vcf_reader):
    test=[]
    for record in vcf_reader:
        if record.num_het==record.num_called:
            continue
        else:
            test.append(record)
        return test