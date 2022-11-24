# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 10:30:14 2021

@author: chris.kerklaan
"""

import sys
import types
import datetime


class Progress:
    def __init__(
        self,
        iterator=None,
        message="",
        quiet=False,
        enum=False,
        total=None,
        yield_progress=False,
        size=30,
        file=sys.stdout,
        start_offset=None,
    ):
        self.iterator = iterator
        self.message = message
        self.size = size
        self.file = file
        self.time = datetime.datetime.now()
        self.enumerate = enum
        self.yield_progress = yield_progress

        if start_offset:
            self.i = start_offset
        else:
            self.i = 0

        if isinstance(iterator, types.GeneratorType) and total == None:
            raise ValueError("Type is generator, please fill in total")

        if total:
            self.total = total
        else:
            self.total = len(iterator)

        self.quiet = quiet

    def __iter__(self):
        if self.enumerate:
            for i, x in enumerate(self.iterator):
                self.show()
                if self.yield_progress:
                    yield self.progress, i, x
                else:
                    yield i, x
        else:
            for i in self.iterator:
                self.show()

                if self.yield_progress:
                    yield self.progress, i
                else:
                    yield i

    def __next__(self):
        yield self.iterator.__next__()

    @property
    def progress(self):
        """progress out of a 100"""
        return self.progress_factor * 100

    @property
    def progress_factor(self):
        return self.i / self.total

    def gdal(self, complete, message=None, cb_data=None):
        """Emit progress report in numbers for 10% intervals and dots for 3%"""
        self.i = int(self.total * complete)
        self.show(False)

    def show(self, add=True, custom_add=0):
        """shows the bar, count can be increased if needed
        params:
            add: automatically adds 1 to the progress
            custom_add: can add more or less than 1 to progress
        """
        if add:
            if custom_add != 0:
                self.i += custom_add
            else:
                self.i += 1

        if self.quiet:
            return

        x = int(self.size * self.progress_factor)

        print(
            f"\r{self.message}: [{'.'*x}{' '*(self.size-x)}] {self.i}/{self.total}",
            end="\r",
            flush=True,
        )

        if self.i >= self.total:
            self.show_end()

    def show_end(self):
        end_time = datetime.datetime.now() - self.time
        # next line
        print("\r")
        # then some message
        print(
            f"\r{self.message}: Completed in {end_time} n = {self.total}",
            end="\r\n",
            flush=True,
        )

    def update(self):
        self.show()
        if self.yield_progress:
            return self.progress
