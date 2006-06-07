#!/usr/bin/env python

# monkeypatch sys.path so we can import Parsing.py
import sys
sys.path.insert(1, '..')
from Parsing import Parser, parsefile

# regular imports now ...
import unittest

# TODO: write the test cases for the known values in input1

if __name__ == '__main__':
    parsefile('input1', True)
