#!/usr/bin/env python

import os
import sys
import re
import copy

filename = sys.argv[1]
with open(filename, 'r') as infile:
    lines = infile.read().splitlines()

match = re.search(r'_(\d)(\d)\.', filename)
N = int(match.group(1))
M = int(match.group(2))

class Contour(object):
    
    class DataPoint(object):
        def __init__(self, text, x, y, z):
            self.text = text
            self.x = x
            self.y = y
            self.z = z
        def __eq__(self, x):
            return self.text == x.text
    
    def __init__(self, label):
        self.label = label
        self.points = []
        
    def add_point(self, *args, **kwargs):
        self.points.append(Contour.DataPoint(*args, **kwargs))


headpat = re.compile(r'^#\s*[Cc]ontour\s+\d+\s*,\s*[Ll]abel\s*:\s*([\-+\.0-9]+)\s*$')
datapat = re.compile(r'^\s*([\-+0-9\.]+)\s+([\-+0-9\.]+)\s+([\-+0-9\.]+)\s*$')
contour = None
contours = []

for line in lines:
    if headpat.match(line):
        label = headpat.match(line).group(1)
        if contour is not None:
            contours.append(contour)
        contour = Contour(label)
    elif datapat.match(line):
        x = datapat.match(line).group(1)
        y = datapat.match(line).group(2)
        z = datapat.match(line).group(3)
        if contour is None:
            contour = Contour(z)
        elif contour.label != z:
            contours.append(contour)
            contour = Contour(z)
        contour.add_point(line, x, y, z)
    elif line == "":
        if contour is not None:
            contours.append(contour)
            contour = None
    elif line.startswith('#'):
        continue
if contour is not None:
    contours.append(contour)

# Filter out any contour objects that have zero points.
contours = [x for x in contours if len(x.points) > 0]


# We have to try glue the contours together into continous lines.
old_contours = contours
contours = []
for co in old_contours:
    merged = False
    for cn in contours:
        if co.points[0] == cn.points[0]:
            cn.points = reversed(co.points) + cn.points
            merged = True
            continue
        elif co.points[0] == cn.points[-1]:
            cn.points = cn.points + co.points
            merged = True
            continue
        elif co.points[-1] == cn.points[0]:
            cn.points = co.points + cn.points
            merged = True
            continue
        elif co.points[-1] == cn.points[-1]:
            cn.points = cn.points + reversed(co.points)
            merged = True
            continue
    if not merged:
        contours.append(co)
            


labelpoints = []
k = 14

low = False
for c in contours:
    if len(c.points) < 2*k+1:
        continue
    if low:
        m = int(float(len(c.points)) * 2./3.)
    else:
        m = int(len(c.points) / 3)
    if c.points[0].z == "1.75":
        if N == 4 and M == 1:
            m = int(float(len(c.points)) * 2.0/5.)
        elif N in [1, 2] and M == 2:
            m = int(float(len(c.points)) * 1.0/8.)
        elif N == 1 and M == 3:
            m = int(float(len(c.points)) * 1.2/8.)
        elif N == 2 and M == 3:
            m = int(float(len(c.points)) * 1.05/8.)
        elif N in [3, 4] and M == 3:
            m = int(float(len(c.points)) * 0.55)
        elif N == 4 and M == 4:
            m = int(float(len(c.points)) * 0.43)
        else:
            m = int(float(len(c.points)) * 1.5/8.)
    if c.points[0].z == "2":
        m = int(float(len(c.points)) * 1.75/8.)
    x = float(c.points[m].x)
    y = float(c.points[m].y)
    if x < 0.03 or 0.775 < x or y < 0.013 or 0.049 < y:
        continue
    low = not low
    labelpoints.append( copy.deepcopy(c.points[m]) )
    for n in xrange(m-k,m+k):
        c.points[n].text = ""
    c.points = c.points[0:m-k] + [c.points[m]] + c.points[m+k:]
    

filename = sys.argv[2]
with open(filename, 'w') as outfile:
    for c in contours:
        outfile.write("\n");
        for p in c.points:
            outfile.write(p.text + "\n");

filename = sys.argv[3]
with open(filename, 'w') as outfile:
    for p in labelpoints:
        outfile.write(p.text + "\n");
