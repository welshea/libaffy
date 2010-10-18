#!/usr/bin/python

import sys

def usage(path):
    print """
%s: delta1[,delta2,...] exprfile1 exprfile2

    Each delta value establishes an additional band of possible expression
    differences such that each delta satisfies one of 0 < X <= delta1,
    delta1 < X <= delta2, and so on.  The size of each band's membership
    is reported in absolute as well as percentage terms.

    exprfile1 and exprfile2 must be of equal dimensions (same number of
    probes, same number of chips).  They should be tab separated format
    as produced by libaffy.
""" % path
    sys.exit()

class bucket:
    def __init__(self, upperbound, relation='<=', isremainder=False):
        self.upperbound = upperbound
        self.count = 0
        self.isremainder = isremainder
        self.relation = relation

def print_buckets(buckets, total):
    for b in buckets:
        print "delta %15s %15u (%04.2f%%)" % (str(b.relation) 
                                             + ' ' +  str(b.upperbound),
                                             b.count,
                                             (1.0*b.count)/total * 100)

def update_buckets(buckets, val):
    for b in buckets:
        if b.isremainder or val <= b.upperbound:
            b.count += 1
            return

if (len(sys.argv) != 4):
    usage(sys.argv[0])

f1 = open(sys.argv[2], mode='r')
f2 = open(sys.argv[3], mode='r')

numdatacols = len(f1.readline().split('\t')) - 1  # first col is file name

if not numdatacols > 0:
    print "ERROR: no data columns or non-tab separated input file"
    sys.exit(1)

if (len(f2.readline().split('\t')) - 1) != numdatacols:
    print "ERROR: column dimensions do not match"
    sys.exit(1)

totalexprs  = 0
bucketlist  = [ bucket(0, relation='=') ]

for i in map(lambda x: float(x), sys.argv[1].split(',')):
    bucketlist += [ bucket(i) ]

bucketlist += [ bucket('', relation='Remainder',  isremainder=True) ]

for line1 in f1:
    line2 = f2.readline()

    if line2 == '':
        print "ERROR: row dimensions do not match"
        sys.exit(1)

    line1vals = map(float, (line1.split('\t')[1:])) # discard probe name
    line2vals = map(float, (line2.split('\t')[1:]))

    for i in xrange(len(line1vals)):
        update_buckets(bucketlist, abs(line1vals[i] - line2vals[i]))

    totalexprs += numdatacols

print "Total number of expression values: " + str(totalexprs) + "\n"

print_buckets(bucketlist, totalexprs)

