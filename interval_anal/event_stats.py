# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
import os


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

class chrom_event_stats_cls(object) :
    file_name_end = '_event_stats.h5'
    table_name_end = '_event_stats'
    def __init__(self, chrom) :
        file_name = 'chrom_' + str(chrom) + self.file_name_end
        table_name = 'chrom_' + str(chrom) + self.table_name_end
        file_path = os.path.join(mod_dir, file_name)
        h5 = tb.open_file(file_path, 'r')
        stats_table = getattr(h5.root, table_name)
        self.stats = stats_table[:]
        h5.close()
