import os
import sys
from skbuild import setup  # This line replaces 'from setuptools import setup'
sys.path.insert(0, os.path.dirname(__file__))
import versioneer


with open("README.md", "r") as fh:
    long_description = fh.read()

configuration = {
    'name': 'sotb-wrapper',
    'packages': ["sotb_wrapper"],
    'package_dir': {'sotb-wrapper': 'sotb_wrapper'},
    'version': versioneer.get_version(),
    'cmdclass': versioneer.get_cmdclass(),
    'description': "wrapper to call fortran routines from SEISCOPE optimization toolbox",
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'url': 'https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes',
    'author': "Oscar Mojica",
    'author_email': 'o_mojical@hotmail.com',
    'include_package_data': True,
    'cmake_args': ["-DSKBUILD=ON"],
    'license': 'MIT',
    'install_requires': ['numpy>=1.20'],
    'zip_safe': False
}

setup(**configuration)
