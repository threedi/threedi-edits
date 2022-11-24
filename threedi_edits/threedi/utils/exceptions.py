# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 09:04:00 2021

@author: chris.kerklaan
"""

_IDK = "¯\_(ツ)_/¯"
_ANGRY = "(╯°□°）╯︵ ┻━┻ ┻━┻ ︵ ヽ(°□°ヽ)"
_HAPPY = "＼(^o^)／"


# warnings and errors
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return f"{_IDK} (Warning):" + str(msg) + "\n"


class ThreediValueTypeError(Exception):
    def __init__(self, message, table, key, value):
        print(_ANGRY)
        self.message = message
        self.table = table
        self.key = key
        self.value = value


class ThreediConnectedTableError(Exception):
    def __init__(self, message, table, connected_table, field, connection_id):
        print(_ANGRY)
        self.message = message
        self.table = table
        self.connected_table = connected_table
        self.connection_id = connection_id
        self.field = field


class ThreediValueTypeWarning(UserWarning):
    def __init__(self, message, table, key, value):
        self.message = message
        self.table = table
        self.key = key
        self.value = value
