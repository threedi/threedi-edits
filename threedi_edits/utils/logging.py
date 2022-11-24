# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:39:17 2019

@author: chris.kerklaan
"""
import logging
import sys


def show_console_logging(level="INFO"):
    if len(logging.root.handlers) == 0:
        if level == "INFO":
            logging.basicConfig(stream=sys.stdout, level=logging.INFO)
        if level == "DEBUG":
            logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
        if level == "WARNING":
            logging.basicConfig(stream=sys.stdout, level=logging.WARNING)
