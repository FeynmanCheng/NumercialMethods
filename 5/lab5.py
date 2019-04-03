# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 09:24:00 2018

@author: Administrator
"""
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import xlrd
if __name__ == "__main__":
    workbook = xlrd.open_workbook(r'C:\lab5.xlsx')
    sheet1 = workbook.sheet_by_index(0)
    data = [[]]
    for i in range(0,4):
        x = sheet1.col_values(2*i)
        y = sheet1.col_values(2*i + 1)
        while(x[-1] == ''):
            x.pop(-1) 
        while(y[-1] == ''):
            y.pop(-1)
        data.append([])
        data[i].append(x)
        data[i].append(y)
    fig = []
    for i in range(0,4):
        p = plt.figure()
        plt.xlabel("Steps")
        plt.ylabel("Accuracy")
        plt.plot(data[i][0],data[i][1])
        plt.xticks(range(0,len(data[i][0])))
        if i < 2:
            title = "Newton method "
        else:
            title = "Secant method "
        if i % 2 == 0:
            title += "x0 = 1"
        else:
            title += "x0 = -2"
        plt.title(title)
        fig.append(p)