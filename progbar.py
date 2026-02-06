#!/usr/bin/env python3

'''
Unicode progress bar for use with "with" statements and "+=" operator.

by Dan Roche, 2015
Should work in Python 2.7 or 3+

Typical usage:

with ProgressBar(total) as pb:
    for x in range(total):
        # do something
        pb += 1

with PBIter([1,2,3,4,5]):
    # do something
'''


from __future__ import print_function, unicode_literals

import math
import subprocess
import sys
import os.path
from operator import length_hint

if sys.version_info < (3,0):
    chr = unichr

class PBIter:
    def __init__(self, collection, size=None, erase=False, basic=False, dest=sys.stderr):
        self._pb = ProgressBar(length_hint(collection) if size is None else size, erase, basic, dest)
        self._iter = iter(collection)
        self._pb.start()

    def __next__(self):
        try:
            res = next(self._iter)
        except StopIteration:
            pass
        else:
            self._pb += 1
            return res
        self._pb.finish()
        raise StopIteration

    def __iter__(self):
        return self

class ProgressBar:
    def __init__(self, maxval=100, erase=False, basic=False, dest=sys.stderr):
        """Creates a new ProgressBar but doesn't display it yet.

        maxval is the value at which the bar will reach 100%
        If erase is True, the bar disappears when it's finished.
        If basic is True, then just write the percentage with no unicode bar characters.
        dest is the file stream to put the bar on.
        """

        self.maxval = maxval
        self.dest = dest
        self.val = 0 # the current value, between 0 and maxval
        self.redraw_at = 0 # past this value the bar will be redrawn

        # partial is an array of unicode characters for partial bars
        self.partial = [' ']
        for code in range(0x258f, 0x2587, -1):
            self.partial.append(chr(code))

        self.gran = len(self.partial) - 1
        self.active = False
        self.basic = basic
        self.erase = erase

    def __del__(self):
        if self.active:
            self.finish()

    def start(self):
        """Starts the bar.
        Typically called from entering a "with" block.
        """
        assert not self.active
        # print(file=self.dest)
        self.active = True
        self.redraw()

    def finish(self):
        """Finishes the bar, even if not at 100%.
        Typically called from leaving a "with" block.
        """
        assert self.active
        self.redraw()
        if self.erase:
            if self.basic:
                nspaces = 14
            else:
                nspaces = self.getwidth()
            print('\r' + ' '*nspaces + '\r', end="", file=self.dest)
        else:
            print('', file=self.dest)
        self.active = False

    def update(self, newval):
        """Changes the stored value and checks if a redraw is necessary.
        Typically called from a += operation.
        """
        assert self.val <= newval <= self.maxval
        self.val = newval
        if self.active and self.val >= self.redraw_at:
            self.redraw()

    def __iadd__(self, increment):
        self.update(self.val + increment)
        return self

    def __int__(self):
        return self.val

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, typ, val, cb):
        self.finish()
        return (typ is None)

    def getwidth(self):
        try:
            termname = os.path.realpath('/dev/stdout')
            return int(subprocess.check_output(['stty', '-F', termname, 'size']).split()[1])
        except subprocess.CalledProcessError:
            # default width if not a terminal
            return 72

    def redraw(self):
        """(re)-draws the bar on the specified output file.
        If not in basic mode, the subprocess is queried to find the terminal width.
        """
        assert self.active
        assert 0 <= self.val <= self.maxval
        percent = self.val * 100 // self.maxval
        if self.basic:
            print('\rProgress: {}%'.format(percent), end="", file=self.dest)
            self.redraw_at = max(self.val+1, ((percent+1)*self.maxval + 99) // 100)
        else:
            total_width = self.getwidth()
            total_blocks = total_width - 6
            total_subblocks = total_blocks * self.gran
            subblocks = self.val * total_subblocks // self.maxval
            nfull, remain = divmod(subblocks, self.gran)
            line = "\r\u2592"
            line += self.partial[-1] * nfull
            if remain:
                line += self.partial[remain]
                line += self.partial[0] * (total_blocks - nfull - 1)
            else:
                line += self.partial[0] * (total_blocks - nfull)
            line += "\u2592{:>3}%".format(percent)
            print(line, end="", file=self.dest)
            self.redraw_at = min(
                ((percent+1)*self.maxval + 99) // 100,
                ((subblocks+1)*self.maxval + total_subblocks-1) // total_subblocks
                )
            assert self.redraw_at > self.val

# test out the different options
def test():
    import random
    import time
    maxval = random.randrange(1234567)
    maxup = 100000
    finish = random.random() < .5
    delay = .1
    with ProgressBar(maxval, random.random() < .5, random.random() < .5) as pbar:
        while int(pbar) < maxval:
            added = random.randrange(maxup)
            if int(pbar) + added > maxval:
                if finish:
                    added = maxval - int(pbar)
                else:
                    break
            pbar += added
            time.sleep(delay)

def usage(error=True):
    print("Usage: {} [OPTIONS] maxval [curval]".format(sys.argv[0]))
    print("  where OPTIONS consists of:")
    print("    -h: Print this help message and exit")
    print("    -r: Do some random testing")
    print("    -e: Erase the bar when finished")
    print("    -b: Make a basic version (no fancy unicode)")
    print("  maxval is the final (finishing) value.")
    print("")
    print("  If curval is given, the state of the bar at that moment")
    print("  is shown and the program exits immediately.")
    print("  Otherwise, the program reads values from stdin.")
    exit(1 if error else 0)

if __name__ == '__main__':
    opts = ""
    maxval = None
    curval = None

    for arg in sys.argv[1:]:
        if arg.startswith('-'):
            opts += arg[1:]
        elif maxval is None:
            try:
                maxval = int(arg)
            except ValueError:
                usage()
        elif curval is None:
            try:
                curval = int(arg)
            except ValueError:
                usage()
        else:
            usage()

    erase = False
    basic = False
    for opt in opts:
        if opt == 'h':
            usage(False)
        elif opt == 'r':
            test()
        elif opt == 'e':
            erase = True
        elif opt == 'b':
            basic = True
        else:
            usage()

    if maxval is None:
        usage()
        exit()

    pb = ProgressBar(maxval, erase, basic)

    if curval is None:
        with pb, sys.stdin as inp:
            for inline in inp:
                newval = int(inline)
                pb.update(newval)
                if newval >= maxval:
                    break
    else:
        pb.update(curval)
        with pb:
            pass
