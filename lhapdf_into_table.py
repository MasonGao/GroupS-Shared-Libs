#!/usr/bin/python
'''
    File: lhapdf_into_table.py
    Author: Kaiyuan Shi
    Date: July 6, 2023
    Description: Generate the table for the given PDF/FF set and save the array in file.
        The tables are given in x, Q in log scale.
    Args:
        name: name of the pdf set, e.g. 'cteq6l1'
        ind: index of set to use, 0 by default
'''

import numpy as np
from numpy import log10, sqrt

from datetime import datetime

import lhapdf
import os
from argparse import ArgumentParser

def PDF_generate(PDF_setname, ind = 0):
    # Get the names of the PDF sets:
    if not os.path.exists('./PDFFF/'):
        os.makedirs('./PDFFF/')

    d_setname = './PDFFF/' + PDF_setname + str(ind) + '_d.txt'
    u_setname = './PDFFF/' + PDF_setname + str(ind) + '_u.txt'
    s_setname = './PDFFF/' + PDF_setname + str(ind) + '_s.txt'
    c_setname = './PDFFF/' + PDF_setname + str(ind) + '_c.txt'
    b_setname = './PDFFF/' + PDF_setname + str(ind) + '_b.txt'
    g_setname = './PDFFF/' + PDF_setname + str(ind) + '_g.txt'

    p = lhapdf.mkPDF(PDF_setname, ind)

    # Get the ranges of x and Q2 for table generation.
    Q2 = np.logspace(log10(p.q2Min), log10(p.q2Max), base = 10, num = 500)
    X = np.logspace(log10(p.xMin), log10(p.xMax), base = 10, num = 500)

    # Begin table generation and record time.
    start = datetime.now()

    print(PDF_setname + str(ind) + ' starting at ' + start.strftime("%Y/%m/%d %H:%M:%S"))
    
    steps_total = len(Q2) * len(X)
    steps_taken = 0
    
    X_arr = np.zeros(steps_total)
    Q_arr = np.zeros(steps_total)
    
    d_arr = np.zeros(steps_total)
    u_arr = np.zeros(steps_total)
    s_arr = np.zeros(steps_total)
    c_arr = np.zeros(steps_total)
    b_arr = np.zeros(steps_total)
    g_arr = np.zeros(steps_total)

    for x in X:
        for q2 in Q2:

            X_arr[steps_taken] = x
            Q_arr[steps_taken] = sqrt(q2)

            d_arr[steps_taken] = p.xfxQ2(1, x, q2)
            u_arr[steps_taken] = p.xfxQ2(2, x, q2)
            s_arr[steps_taken] = p.xfxQ2(3, x, q2)
            c_arr[steps_taken] = p.xfxQ2(4, x, q2)
            b_arr[steps_taken] = p.xfxQ2(5, x, q2)
            g_arr[steps_taken] = p.xfxQ2(21, x, q2)

            steps_taken += 1


    d_table = np.zeros((len(X_arr), 3))
    u_table = np.zeros((len(X_arr), 3))
    s_table = np.zeros((len(X_arr), 3))
    c_table = np.zeros((len(X_arr), 3))
    b_table = np.zeros((len(X_arr), 3))
    g_table = np.zeros((len(X_arr), 3))

    for i in range(len(X_arr)):
        d_table[i] = [X_arr[i], Q_arr[i], d_arr[i]]
        u_table[i] = [X_arr[i], Q_arr[i], u_arr[i]]
        s_table[i] = [X_arr[i], Q_arr[i], s_arr[i]]
        c_table[i] = [X_arr[i], Q_arr[i], c_arr[i]]
        b_table[i] = [X_arr[i], Q_arr[i], b_arr[i]]
        g_table[i] = [X_arr[i], Q_arr[i], g_arr[i]]

    # Also generate the alpha_s table, which is an 1D array.
    aQ2 = np.logspace(log10(0.5), 11, base = 10, num = 2000)
    a_arr = np.zeros(len(aQ2))
    a_setname = '../PDFFF/' + PDF_setname + str(ind) + '_a.txt'

    steps_taken = 0
    for aq2 in aQ2:
        a_arr[steps_taken] = p.alphasQ2(aq2)
        steps_taken += 1

    a_table = np.zeros((len(a_arr), 2))
    for i in range(len(a_arr)):
        a_table[i] = [aQ2[i], a_arr[i]]

    # Save the tables.
    np.savetxt(d_setname, d_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(u_setname, u_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(s_setname, s_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(c_setname, c_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(b_setname, b_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(g_setname, g_table, delimiter=' & ', header = 'x & Q & F', newline='\n', fmt='%.16f')
    np.savetxt(a_setname, a_table, delimiter=' & ', header = 'Q2 & a', newline='\n', fmt='%.16f')

    end = datetime.now()
    print(PDF_setname + str(ind) + ' completed at ' + end.strftime("%Y/%m/%d %H:%M:%S") + ', using ' + str(end - start))

parser = ArgumentParser(description = 'Generate the PDF or FF datasets given the name and index.')
parser.add_argument('--name', type = str, help = 'The name of the dataset for table generation.')
parser.add_argument('--index', type = int, default = 0, help = 'The index of dataset')
args = parser.parse_args()

if __name__ == "__main__":

    PDF_setname = args.name
    ind = int(args.index)
    PDF_generate(PDF_setname, ind)
