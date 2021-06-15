"""Path hack to make tests work."""

import os
import sys

bp = os.path.dirname(os.path.realpath('.')).split(os.sep)
modpath = os.sep.join(bp[:-1] + ['seiscope_opt_tb_wrapper'])
sys.path.insert(0, modpath)
