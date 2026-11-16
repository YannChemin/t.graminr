"""Put the addon root on sys.path so the tests import graminr_lib directly.

SPDX-License-Identifier: Unlicense
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
