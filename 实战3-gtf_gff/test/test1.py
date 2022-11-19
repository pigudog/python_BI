# -*- coding: utf-8 -*-
"""
日期：19/11/2022 13:33
作者： PigUDog
"""
import os
import sys
from os import path
d = path.dirname(__file__)
parent_path = os.path.dirname(d)
sys.path.append(parent_path)
import gtftools # sys.path
def test_info():
    # 测试一些基本信息
    print(gtftools.__version__)
    print(gtftools.__author__)
    print(gtftools.__email__)

def test_gff3():
    gff3 = Gff3(file="../data/test.gff3")
    gff3.to_gtf(file="./test_to_gtf.gtf")


def main():
    # test_info()
    pass

if __name__ == '__main__':
    main()
