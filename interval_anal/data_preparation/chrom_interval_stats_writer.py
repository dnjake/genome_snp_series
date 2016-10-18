# -*- coding: utf-8 -*-

import tables as tb
from interval_stats import chrom_interval_stats_cls

class chrom_stats_writer_cls(object) :
    file_name_end = '_interval_stats.h5'
    table_name_end = '_interval_stats'
        
    def do_chrom(self, chrom) :
        print 'writing stats for chrom', chrom                
        chrom_stats = chrom_interval_stats_cls(chrom)
        chrom_stats.process_intervals()
        stats = chrom_stats.stats
        out_file_name = 'chrom_' + str(chrom) + self.file_name_end
        out_table_name = 'chrom_' + str(chrom) + self.table_name_end
        h5 = tb.open_file(out_file_name, 'w')
        out_data = h5.create_table('/', out_table_name, description=stats.dtype)
        out_data.append(stats)
        h5.close()
        
csw = chrom_stats_writer_cls()
for chrom in range(1,23) :
    csw.do_chrom(chrom)
    
    