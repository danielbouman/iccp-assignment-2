"""
This function returns a message to the user in the console. This message is intended to notify the user that the simulation has
started running, and is most likely still running. This function has been optimized for speed on quadcore computers.
"""
from __future__ import print_function       # make print function work in older versions of python
def message():
    print(" ")
    print(" ")
    print("* _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ *")
    print("~_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _~")
    print(" _______   ___  ___ ___   __ ___   __ ____  ___   __ ________                  ")
    print("/|    _ | ||  | |  ||  |  | ||  |  | ||   | |  |  | ||       |                 ")
    print("||   | || ||  | |  ||   |_| ||   |_| ||   | |   |_| ||    ___|                 ")
    print("||   |_||_||  |_|  ||       ||       ||   | |       ||   |___                  ")
    print("||    __  ||       ||  _    ||  _    ||   | |  _    ||   ||  |____  ____  ____ ")
    print("||   | || ||       || | |   || | |   ||   | | | |   ||   |_| ||   | |   | |   |")
    print("||___| ||_||_______||_|  |__||_|  |__||___| |_|  |__||_______||___| |___| |___|")
    print("|___/  |_/|_______/|_/   |_/ |_/  ||   | | |  |  ||       ||   | |   | |   |")
    print("~_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _~")
    print("* _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ *")
    print(" ")
    print(" ")
    return