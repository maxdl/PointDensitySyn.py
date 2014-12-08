from distutils.core import setup
import py2exe
import version

setup(
    windows = [
        {
            "script": "PointDensitySyn.py",
            "icon_resources": [(1, version.title)]
        }
    ],
)

