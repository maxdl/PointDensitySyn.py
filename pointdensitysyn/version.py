# -*- coding: utf-8 -*-

import os.path
import sys

title = "PointDensitySyn"
author = "Max Larsson"
version = "1.0"
date = ("October", "6", "2016")
email = "max.larsson@liu.se"
homepage = "www.hu.liu.se/forskning/larsson-max/software"
if hasattr(sys, 'frozen'):
    if '_MEIPASS2' in os.environ: 
        path = os.environ['_MEIPASS2']
    else:
        path = sys.argv[0]
else:
    path = __file__   
app_path = os.path.dirname(path)
icon = os.path.join(app_path, "pds.ico")

