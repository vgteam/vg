#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
blockify.py: make sure input and output FDs are in blocking mode before running a command.
"""

import os
import sys
import fcntl

# Check input
if len(sys.argv) < 2:
    raise RuntimeError("A program to execute is required")

# Fix up standard file descriptors (0, 1, 2) by clearing the nonblocking flag.
# See https://stackoverflow.com/a/30172682
for fd in range(3):
    fcntl.fcntl(fd, fcntl.F_SETFL, fcntl.fcntl(fd, fcntl.F_GETFL) & ~os.O_NONBLOCK)

# Become the first argument running with the rest of the arguments.
# Make sure to pass along the program name.
os.execvp(sys.argv[1], sys.argv[1:])
