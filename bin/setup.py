#!/usr/bin/env python

import os
import sys
import path
import json

from libcommon import *

HOST_HOME = path.Dir(HOST_HOME, build=True)
with JOBs_json.open("wt") as fout:
    fout.write(json.dumps([], indent=2))
with HOSTs_json.open("wt") as fout:
    fout.write(json.dumps([], indent=2))
