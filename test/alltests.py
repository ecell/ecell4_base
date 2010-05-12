#!/usr/bin/env python

import unittest
import os.path
import glob

def suite():

    test_files = glob.glob('*_test.py')
    modules_to_test = [os.path.splitext(file)[0] for file in test_files]

    alltests = unittest.TestSuite()

    for module in map(__import__, modules_to_test):
        alltests.addTest(unittest.findTestCases(module))

    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
