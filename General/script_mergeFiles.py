# -*- coding: latin-1 -*-.

#! /usr/bin/env python

from __future__ import with_statement

f1 = open('FEDId_FEDCh_APVId_AverageAPVCMN-Minus127_PbPb2024_Run387799_2sigma.txt')
f2 = open('FEDId_FEDCh_APVId_APVsWithTwoSigmasOutInMisbehaving_AllFEDs.txt')
f3 = open('output.txt','w')

column_separator = ' '  # one space to separate the columns

array_text = []
for line1 in f1.readlines():
    columns1 = line1.split(column_separator)  # Get the data groups
    f2.seek(0)
    for line2 in f2.readlines():
        columns2 = line2.split(column_separator)  # Get the data groups
        if columns1[0]==columns2[0] and columns1[1]==columns2[1] and columns1[2]==columns2[2]:
            line1 = line1.replace(columns1[3],columns2[3])
    array_text.append(line1)
f3.writelines(array_text)
f1.close()
f2.close()
f3.close()
