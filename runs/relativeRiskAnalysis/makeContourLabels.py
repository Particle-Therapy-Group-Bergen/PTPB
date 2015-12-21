#!/usr/bin/env python

import os
import sys
import re
import copy

filename = sys.argv[1]
with open(filename, 'r') as infile:
    lines = infile.read().splitlines()

match = re.search(r'contours_(\w+)_(\d)(\d)\.', filename)
organtech = match.group(1)
N = int(match.group(2))
M = int(match.group(3))

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
k = 16

low = False
for c in contours:
    if len(c.points) < 2*k+1:
        continue
    if low:
        m = int(float(len(c.points)) * 2./3.)
    else:
        m = int(len(c.points) / 3)
    # The following is for adjusting individual contour labels.
    # N indicates the plot column left to right 1 .. 4.
    # M indicates the plot row top to bottom 1 .. 4.
    # Here we adjust the distance along the curve m (in the range 0..1) to place
    # the label.
    if organtech == "cion_bladder":
        if c.points[0].z == "1.75":
            if N == 4 and M == 1:
                m = int(float(len(c.points)) * 2.0/5.)
            elif N == 1 and M == 2:
                m = int(float(len(c.points)) * 0.14)
            elif N == 2 and M == 2:
                m = int(float(len(c.points)) * 0.115)
            elif N == 3 and M == 2:
                m = int(float(len(c.points)) * 0.19)
            elif N == 4 and M == 2:
                m = int(float(len(c.points)) * 0.19)
            elif N == 1 and M == 3:
                m = int(float(len(c.points)) * 0.21)
            elif N == 2 and M == 3:
                m = int(float(len(c.points)) * 0.17)
            elif N in [3, 4] and M == 3:
                m = int(float(len(c.points)) * 0.31)
            elif N == 1 and M == 4:
                m = int(float(len(c.points)) * 0.2)
            elif N == 2 and M == 4:
                m = int(float(len(c.points)) * 0.185)
            elif N == 3 and M == 4:
                m = int(float(len(c.points)) * 0.22)
            elif N == 4 and M == 4:
                m = int(float(len(c.points)) * 0.43)
            else:
                m = int(float(len(c.points)) * 1.5/8.)
        if c.points[0].z == "1.5":
            if M == 3:
                m = int(float(len(c.points)) * 1./2.)
        if c.points[0].z == "2":
            if M == 4:
                if N == 1:
                    m = int(float(len(c.points)) * 0.165)
                elif N == 4:
                    m = int(float(len(c.points)) * 0.14)
                else:
                    m = int(float(len(c.points)) * 0.15)
            else:
                if N in [1, 2]:
                    m = int(float(len(c.points)) * 1.75/8.)
                elif N == 3:
                    m = int(float(len(c.points)) * 1.65/8.)
                else:
                    m = int(float(len(c.points)) * 1.5/8.)
    elif organtech == "cion_rectum":
        if float(c.points[0].z) > 1:
            m = int(float(len(c.points)) * 0.0)
    elif organtech == "proton_bladder":
        if c.points[0].z == "2.5":
            if M == 1 and N == 2:
                m = int(float(len(c.points)) * 0.8)
            elif M == 1 and N == 3:
                m = int(float(len(c.points)) * 0.63)
            elif M == 1 and N == 4:
                m = int(float(len(c.points)) * 0.57)
            elif M == 2 and N == 3:
                m = int(float(len(c.points)) * 0.6)
            elif M == 2 and N == 4:
                m = int(float(len(c.points)) * 0.43)
            elif M == 3 and N == 2:
                m = int(float(len(c.points)) * 0.0)
            elif M == 3 and N == 3:
                m = int(float(len(c.points)) * 0.43)
            elif M == 3 and N == 4:
                m = int(float(len(c.points)) * 0.32)
            elif M == 4 and N in [3, 4]:
                m = int(float(len(c.points)) * 0.32)
    elif organtech == "proton_rectum":
        if float(c.points[0].z) > 1.5:
            m = int(float(len(c.points)) * 0.0)
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
        # Change label for the "1" contour to "1.0".
        if '.' not in p.z and int(p.z) == 1:
            p.text = "{0}  {1}  {2:.1f}".format(p.x, p.y, float(p.z))
        # Ajust the middle position of some of the labels:
        if organtech == "cion_bladder":
            if p.z == "1.75" and M == 2 and N == 1:
                p.text = "{0}  {1}  {2}".format(float(p.x) - 0.005, p.y, p.z)
            elif p.z == "1.75" and M == 2 and N == 2:
                p.text = "{0}  {1}  {2}".format(float(p.x) - 0.0075, p.y, p.z)
            elif p.z == "1.75" and M == 2 and N in [3, 4]:
                p.text = "{0}  {1}  {2}".format(float(p.x) - 0.01, p.y, p.z)
        elif organtech == "cion_rectum":
            if p.z == "0.75" and M == 2:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
            if p.z == "0.75" and M in [3, 4]:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.015, p.y, p.z)
            elif p.z == "1":
                p.text = "{0}  {1}  {2:.1f}".format(float(p.x) + 0.01, p.y, float(p.z))
        elif organtech == "proton_bladder":
            if p.z == "2.5" and M == 1 and N == 2:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.018, p.y, p.z)
            elif p.z == "2.5" and M == 1 and N == 3:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.012, p.y, p.z)
            elif p.z == "2.5" and M == 1 and N == 4:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
            elif p.z == "2.5" and M == 2 and N == 3:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.014, p.y, p.z)
            elif p.z == "2.5" and M == 2 and N == 4:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
            elif p.z == "2.5" and M == 3 and N == 3:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.015, p.y, p.z)
            elif p.z == "2.5" and M == 3 and N == 4:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.008, p.y, p.z)
            elif p.z == "2.5" and M == 4 and N in [1]:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.024, p.y, p.z)
            elif p.z == "2.5" and M == 4 and N in [2, 3]:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.015, p.y, p.z)
            elif p.z == "2.25" and M == 1:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
            elif p.z == "2.25" and M == 2 and N == 4:
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
        elif organtech == "proton_rectum":
            if p.z == "1.5":
                p.text = "{0}  {1}  {2}".format(float(p.x) + 0.01, p.y, p.z)
        outfile.write(p.text + "\n");
