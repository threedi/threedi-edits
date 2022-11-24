# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:56:20 2019

@author: chris.kerklaan - N&S
"""
# First-party imports
import os
import sys
import json
import logging
import datetime
import shutil
import tempfile
from pathlib import Path

from inspect import getmembers, isfunction, isclass, ismodule

# Globals
time_last_print_small = datetime.datetime.now()
time_last_print_large = datetime.datetime.now()


def log_time(log_type, s1="", s2="", size="s"):
    print_time(s1=s1, s2=s2, size=size)

    if log_type == "info":
        logging.info("%s - %s", s1, s2)
    elif log_type == "debug":
        logging.debug("%s - %s", s1, s2)
    elif log_type == "warning":
        logging.warning("%s - %s", s1, s2)
    elif log_type == "error":
        logging.error("%s - %s", s1, s2)
    else:
        print("no log")


class Logger(object):
    def __init__(self, path, quiet=False):
        # self.terminal = sys.stdout
        self.log = open(path, "a")
        self.quiet = quiet
        self.path = path

    def set_subject(self, subject):
        self.s1 = subject

    def show(self, s2="", quiet=None, size="s"):
        # if self.quiet or quiet:
        #     be_quiet=True
        # else:
        #     be_quiet=False

        # if not be_quiet:
        global time_last_print_small
        global time_last_print_large

        now = datetime.datetime.now()
        if size == "s":
            _print = "{} - {} - {} - {}".format(
                now, now - time_last_print_small, self.s1, s2
            )

        elif size == "l":
            _print = "{} - {} - {} - {}".format(
                now, now - time_last_print_large, self.s1, s2
            )

        else:
            pass

        time_last_print_small = now
        time_last_print_large = now

        print("\n" + _print + "\n")
        self.log.write(_print + "\n")

    def close(self):
        self.log = None

    def write(self):
        self.log = None
        self.log = open(self.path, "a")


def print_time(s1="", s2="", size="s", quiet=False):
    if not quiet:
        global time_last_print_small
        global time_last_print_large

        now = datetime.datetime.now()
        if size == "s":
            _print = "{} - {} - {} - {}".format(
                now, now - time_last_print_small, s1, s2
            )

        elif size == "l":
            _print = "{} - {} - {} - {}".format(
                now, now - time_last_print_large, s1, s2
            )

        else:
            pass

        time_last_print_small = now
        time_last_print_large = now
        return print(_print + "\n")


def percentage(count, total):
    return str((count / total) * 100) + "%"


def mk_temp(path=os.getcwd()):
    tempfolder = os.path.join(path, "temp")
    if not os.path.exists(tempfolder):
        os.mkdir(tempfolder)
    return tempfolder


def print_list(_list, subject):
    if not isinstance(_list, list):
        _list = [_list]

    print("\n {}:".format(subject))
    for path in _list:
        print("\t{}".format(str(path)))


def print_dictionary(_dict, subject):
    print("\n {}:".format(subject))
    for key, value in _dict.items():
        print("\t{}:\t\t{}".format(str(key), str(value)))


def write_dictionary(path, _dict):
    with open(path, "w") as write_file:
        json.dump(_dict, write_file, indent=4)


def load_dictionary(path):
    with open(path) as json_file:
        data = json.load(json_file)
    return data


def mk_dir(path=os.getcwd(), folder_name="temp", overwrite=True):
    tempfolder = os.path.join(path, folder_name)

    if os.path.exists(tempfolder) and overwrite:
        shutil.rmtree(tempfolder)
        os.mkdir(tempfolder)
    elif not os.path.exists(tempfolder):
        os.mkdir(tempfolder)
    else:
        pass  # Folder exits, path
    return tempfolder


def create_temporary_copy(path):
    temp = tempfile.NamedTemporaryFile("w+t", delete=True)
    shutil.copy2(temp.name, path)
    return temp.name


class Progress:
    def __init__(self, total, message, size=30, file=sys.stdout):
        self.total = total
        self.message = message
        self.size = size
        self.file = file
        self.i = 0
        self.time = datetime.datetime.now()

        # clear the current line
        print("\r")

    def update(self, quiet=False):
        if quiet:
            return

        self.i += 1
        x = int(self.size * self.i / self.total)

        print(
            f"\r{self.message}: [{'.'*x}{' '*(self.size-x)}] {self.i}/{self.total}",
            end="\r",
            flush=True,
        )

        if self.i >= self.total:
            end_time = datetime.datetime.now() - self.time
            print(
                f"\r{self.message}: Completed {self.total} loops in {end_time}",
                end="\r",
                flush=True,
            )


def pretty_dict(d, indent=0):
    string = ""
    for key, value in d.items():
        string = string + "\t" * indent + str(key)
        if isinstance(value, dict):
            pretty_dict(value, indent + 1)
        else:
            string = string + "\t" * (indent + 1) + str(value)
    return string


class DisplayablePath(object):
    display_filename_prefix_middle = "├──"
    display_filename_prefix_last = "└──"
    display_parent_prefix_middle = "    "
    display_parent_prefix_last = "│   "

    def __init__(self, path, parent_path, is_last):
        self.path = Path(str(path))
        self.parent = parent_path
        self.is_last = is_last
        if self.parent:
            self.depth = self.parent.depth + 1
        else:
            self.depth = 0

    @property
    def displayname(self):
        if self.path.is_dir():
            return self.path.name + "/"
        return self.path.name

    @classmethod
    def make_tree(cls, root, parent=None, is_last=False, criteria=None):
        root = Path(str(root))
        criteria = criteria or cls._default_criteria

        displayable_root = cls(root, parent, is_last)
        yield displayable_root

        children = sorted(
            list(path for path in root.iterdir() if criteria(path)),
            key=lambda s: str(s).lower(),
        )
        count = 1
        for path in children:
            is_last = count == len(children)
            if path.is_dir():
                yield from cls.make_tree(
                    path, parent=displayable_root, is_last=is_last, criteria=criteria
                )
            else:
                yield cls(path, displayable_root, is_last)
            count += 1

    @classmethod
    def _default_criteria(cls, path):
        return True

    def displayable(self):
        if self.parent is None:
            return self.displayname

        _filename_prefix = (
            self.display_filename_prefix_last
            if self.is_last
            else self.display_filename_prefix_middle
        )

        parts = ["{!s} {!s}".format(_filename_prefix, self.displayname)]

        parent = self.parent
        while parent and parent.parent is not None:
            parts.append(
                self.display_parent_prefix_middle
                if parent.is_last
                else self.display_parent_prefix_last
            )
            parent = parent.parent

        return "".join(reversed(parts))


class Files:
    def __init__(
        self,
        file=None,
        root=None,
        folder_only=False,
        folder_filters=["__pycache__", "data"],
        file_filters=[".pyc"],
    ):

        if file:
            root = Path(file).parents[0]

        self.paths = DisplayablePath.make_tree(root)
        self.filter(folder_only, file_filters, folder_filters)

    def __repr__(self):
        return self.general

    def filter(self, folder_only, files=[], folders=[]):

        self.general = ""
        for path in self.paths:
            path_display = path.displayable()

            # flters on folders only
            if folder_only:
                if "." in path_display:
                    continue

            # filters on files
            do_continue = False
            for file in files:
                if path_display.endswith(file):
                    do_continue = True

            if do_continue:
                continue

            do_continue = False

            parents = []
            parent = path.parent

            # if one of the upper parents is that folder, it is skipped
            while parent != None:
                if parent:
                    parents.append(parent.displayname)

                parent = parent.parent

            do_continue = False
            for folder in folders:
                do_continue = len([p for p in parents if p == folder + "/"]) > 0

                if do_continue:
                    break

            for folder in folders:
                if folder + "/" == path.displayname:
                    do_continue = True
                    break

            if do_continue:
                continue

            self.general = self.general + path.displayable() + "\n"


class Functions:
    """retrieves all functions of a script"""

    def __init__(self, module_name):
        self.module = sys.modules[module_name]

    def __repr__(self):
        members = "Functions: \n"
        for name, a in getmembers(self.module, isfunction):
            members = members + "\t-\t" + name + "\n"

        return members


class Classes:
    """retrieves all classes"""

    def __init__(self, module_name, local_only=False):
        self.module_name = module_name
        self.module = sys.modules[module_name]
        self.local_only = local_only

    def __repr__(self):
        members = "Classes: \n"
        for name, _cls in getmembers(self.module, isclass):
            if self.local_only:
                if _cls.__module__ == self.module_name:
                    members = members + "\t-\t" + name + "\n"
            else:
                members = members + "\t-\t" + name + "\n"

        return members


class Modules:
    """retrieves all classes"""

    def __init__(self, module_name):
        self.module = sys.modules[module_name]

    def __repr__(self):
        members = "Modules: \n"
        for name, a in getmembers(self.module, ismodule):
            if not name.isupper():
                members = members + "\t-\t" + name + "\n"

        return members
