# -*- coding: utf-8 -*-

import numpy as np
from autosome_snp_data.chrom_snp_series_data_rdr import chrom_snp_series_data_tables_cls


class chrom_interval_zero_pos_cls(object) :
    event_dtype = np.dtype([('event_pos', 'u4'), ('is_first_pos', '?')])
    interval_dtype = np.dtype([('interval_first_pos', 'u4'), ('interval_last_pos', 'u4')])
    min_snp_cond = cond = 'item_count > 3'
    def __init__(self, chrom) :
        self.chrom = chrom
        self.read_chrom_series_data()
        self.find_intervals()
    
    def read_chrom_series_data(self) :
        rdr = chrom_snp_series_data_tables_cls(self.chrom)
        sdt = rdr.series_data_table
        data_indexes = sdt.read_where(self.min_snp_cond, field='data_index')
        data_size = data_indexes.size
        events = np.zeros(2*data_size, self.event_dtype )
        events['event_pos'][:data_size] = sdt.read_coordinates(data_indexes, field='first_pos')
        events['is_first_pos'][:data_size] = True
        events['event_pos'][data_size:] = sdt.read_coordinates(data_indexes, field='last_pos')
        rdr.close()
        self.events = events

    def find_intervals(self) :
        self.events.sort(order='event_pos')        
        series_count = 0
        intervals = []
        first_pos = None
        for event in self.events :
            if event['is_first_pos'] :
                series_count += 1
                if first_pos is None :
                    first_pos = event['event_pos']
            else :
                series_count -= 1
                if series_count == 0 :
                    last_pos = event['event_pos']
                    intervals.append((first_pos, last_pos))
                    first_pos = None                
        self.series_intervals = np.array(intervals, self.interval_dtype)






