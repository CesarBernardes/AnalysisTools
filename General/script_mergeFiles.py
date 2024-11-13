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


###to add missing APVs
f2.seek(0)
for line3 in f2.readlines():
    foundAPV = False
    columns3 = line3.split(column_separator)  # Get the data groups
    f1.seek(0)
    for line4 in f1.readlines():
        columns4 = line4.split(column_separator)  # Get the data groups
        if columns3[0]==columns4[0] and columns3[1]==columns4[1] and columns3[2]==columns4[2]:
            foundAPV = True
    if foundAPV==False:
        array_text.append(line3)

f3.writelines(array_text)
f1.close()
f2.close()
f3.close()
