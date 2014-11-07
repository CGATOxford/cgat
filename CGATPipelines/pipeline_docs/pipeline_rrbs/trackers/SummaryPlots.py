import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from rrbsReport import *


class readPositionMethylationBias(TrackerImages):
    glob = "plots.dir/*_read_position_methylation_bias.png"
