"""Path hack to make tests work."""

import os
import sys

bp = os.path.dirname(os.path.realpath('.')).split(os.sep)
modpath = os.sep.join(bp + ['seiscope_opt_toolbox_w_ctypes/seiscope_opt_tb_wrapper'])
print(bp, modpath)
sys.path.insert(0, modpath)
