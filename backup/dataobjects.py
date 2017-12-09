#
# -*- coding: utf-8 -*-
# GitHub  dataobjects
# =============================================================================
# Name                : dataobjects.py
# Date                : 2017-08-29 23:31:59
# Author              : sol courtney
# GitHub              : https://github.com/swc2124
# Affiliation         : Columbia University NYC, NY
# Email               : swc2124@columbia.edu
# Language            : Python
# Last Modified by    : swc21
# Last Modified time  : 2017-08-30 00:18:15
# =============================================================================


class DataRow(list):
    """Simple object to serve as the list of data items that go into a
    single row of the output table for skysearcher.py"""
    def __init__(self, table):
        super(DataRow, self).__init__()
        self.table = table
