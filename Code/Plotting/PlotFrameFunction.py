#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import os



# cmdargs = list(map(str,sys.argv))

# filename =  cmdargs[1]

# outputfilename = os.path.splitext(filename)[0]+'.jpg'

# fig = plt.figure()

# ax = fig.add_subplot(111,projection='3d')

def graph(filename):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    outputfilename = os.path.splitext(filename)[0]+'.jpg'
    with open(filename, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        for eachLine in spamreader:
            xpos = eachLine[0]
            ypos = eachLine[1]
            zpos = eachLine[2]

            ax.scatter(xpos,ypos,zpos,color='blue')

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

    plt.savefig(outputfilename, format='jpg')
    plt.close(fig)


def graph(filename,colarray):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    outputfilename = os.path.splitext(filename)[0]+'.png'
    with open(filename, newline='') as csvfile:
        it = 0
        spamreader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        for eachLine in spamreader:
            xpos = eachLine[0]
            ypos = eachLine[1]
            zpos = eachLine[2]
            if(it>=len(colarray)) :
                it = 0

            ax.scatter(xpos,ypos,zpos,c=colarray[it],s=1)
            it = it + 1

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

    plt.savefig(outputfilename, format='png', dpi=1200)
    plt.close(fig)


def graph(filename,colarray,sa):
    fig = plt.figure( figsize=(4, 3) )
    ax = fig.add_subplot(111,projection='3d')
    outputfilename = os.path.splitext(filename)[0]+'.png'
    with open(filename, newline='') as csvfile:
        it = 0
        spamreader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        for eachLine in spamreader:
            xpos = eachLine[0]
            ypos = eachLine[1]
            zpos = eachLine[2]
            if(it>=len(colarray)) :
                it = 0
            ss = 1
            if(it<sa) :
                ss = 10
            ax.scatter(xpos,ypos,zpos,c=colarray[it],s=ss)
            it = it + 1

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

    plt.savefig(outputfilename, format='png', dpi=2400)
    plt.close(fig)