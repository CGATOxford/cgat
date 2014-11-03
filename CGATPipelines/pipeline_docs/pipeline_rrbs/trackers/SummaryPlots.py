import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from rrbsReport import *


class TrackerImages(Trackerimages):
    glob = "plots.dir/*_read_position_methylation_bias.png"

