#!/usr/bin/python
from nbody6tools import __utils
import sys

parser = __utils.executable_options()
args = parser.parse_args()

__utils.execute_actions(args)
