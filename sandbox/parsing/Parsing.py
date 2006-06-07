#!/usr/bin/env python
"""
This module provides a similar interface to Citcom's input-file parser.
"""

# TODO: use properly named exception subclasses

def strip_comments(line):
    """Strip all shell-style comments (first '#' to EOL) from line."""
    k = line.find('#')
    if k >= 0:
        return line[:k]
    return line

def strip_quotes(s):
    """Given a quoted string, strip the quotes."""
    if s and ((s[0] == '"' == s[-1]) or (s[0] == "'" == s[-1])):
        return s[1:-1]
    return s

def parse_eqn(stmt):
    """Given an assignment statement, return the
    parameter name and parameter value."""
    k = stmt.find('=')
    if k < 0:
        # XXX: use a properly named exception
        raise Exception('%r is not an assignment!' % stmt)
    name, value = stmt[:k], stmt[(k+1):]
    return (name, value)


# simple function to return a dictionary
def parsefile(filename, verbose=False):
    """Returns dictionary of parameter names, whose values
    are taken from the specified input file."""
    import sys
    params = dict()
    contents = open(filename, 'r').readlines()
    if verbose:
        print >> sys.stderr, '# Parsing input-file %r' % filename
    for line in contents:
        line = strip_comments(line.strip())
        if not line:
            continue
        assignments = line.split()
        for assign in assignments:
            name, value = parse_eqn(assign)
            if verbose:
                print >> sys.stderr, 'params[%r] = %r' % (name, value)
            params[name] = strip_quotes(value)
    return params


def input_method(cast, vectorize=False):
    """Generate an appropriate input method for the Parser wrapper class."""
    def input_wrapper(self, name, **kw):        
        if name not in self.params:
            if kw.get("essential", False):
                raise Exception("Required parameter %r is missing!" % name)
            return kw.get("default", None)
        if vector:
            values = self.params[name].split(',')
            return [ cast(x) for x in values ]
        value = cast(self.params[name])
        if cast in (int, float):
            minval = kw.get("minvalue")
            maxval = kw.get("maxvalue")
            if minval and value:
                assert float(minval) <= value
            if maxval and value:
                assert float(maxval) >= value
        return value
    return input_wrapper

# finally, a wrapper object
class Parser(object):
    """Wrapper class around the parsefile() dictionary.
    It emulates the interface of Citcom's input file parser."""

    def __init__(self, filename, verbose=False):
        self.params = parsefile(filename, verbose)
    
    input_str    = input_method(str)
    input_bool   = input_method(bool)
    input_int    = input_method(int)
    input_float  = input_method(float)
    input_double = input_method(float)

    input_str_vector    = input_method(str, vectorize=True)
    input_bool_vector   = input_method(bool, vectorize=True)
    input_int_vector    = input_method(int, vectorize=True)
    input_float_vector  = input_method(float, vectorize=True)
    input_double_vector = input_method(float, vectorize=True)


if __name__ == '__main__':
    parsefile('test/input1', True)
