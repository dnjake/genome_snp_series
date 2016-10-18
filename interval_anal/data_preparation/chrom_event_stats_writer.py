# -*- coding: utf-8 -*-

import tables as tb
from chrom_event_stats import chrom_event_stats_cls

class chrom_events_writer_cls(object) :
    file_name_end = '_event_stats.h5'
    table_name_end = '_event_stats'
        
    def do_chrom(self, chrom) :
        print 'writing events for chrom', chrom                
        event_stats = chrom_event_stats_cls(chrom)
        event_stats.process_intervals()
        stats = event_stats.stats
        out_file_name = 'chrom_' + str(chrom) + self.file_name_end
        out_table_name = 'chrom_' + str(chrom) + self.table_name_end
        h5 = tb.open_file(out_file_name, 'w')
        out_data = h5.create_table('/', out_table_name, description=stats.dtype)
        out_data.append(stats)
        h5.close()
        
cew = chrom_events_writer_cls()
for chrom in range(1,23) :
    cew.do_chrom(chrom)
