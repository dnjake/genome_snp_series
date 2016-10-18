# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
import os


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

class chrom_interval_stats_cls(object) :
    ratio_dtype = np.dtype([('index', 'i4'), ('ratio', 'f4')])
    zero_interval_dtype = np.dtype([('end_last', 'u4'), ('start_next', 'u4'), ('zero_length', 'u4')])
    file_name_end = '_interval_stats.h5'
    table_name_end = '_interval_stats'
    def __init__(self, chrom) :
        file_name = 'chrom_' + str(chrom) + self.file_name_end
        table_name = 'chrom_' + str(chrom) + self.table_name_end
        file_path = os.path.join(mod_dir, file_name)
        h5 = tb.open_file(file_path, 'r')
        stats_table = getattr(h5.root, table_name)
        self.stats = stats_table[:]
        h5.close()
        
    def array_for_stats_ratio(self) :
        data = np.zeros(self.stats.size, self.ratio_dtype)
        data['index'] = np.arange(self.stats.size)
        return data
        
    def stats_by_grouped_snps_to_interval_length(self, return_ratio=False)  :
        ratio_data = self.array_for_stats_ratio()
        grouped_snps = self.stats['grouped_snp_count'].astype('f4')
        interval_lengths = self.stats['interval_length'].astype('f4')
        ratio_data['ratio'] = grouped_snps/interval_lengths
        ratio_data.sort(order='ratio')
        ratio_data = ratio_data[::-1]
        out_stats = self.stats[ratio_data['index']]
        if return_ratio :
            return out_stats, ratio_data
        else :
            return out_stats

    def stats_series_to_interval_length(self, return_ratio=False)  :
        ratio_data = self.array_for_stats_ratio()
        series_counts = self.stats['interval_series_count'].astype('f4')
        interval_lengths = self.stats['interval_length'].astype('f4')
        ratio_data['ratio'] = series_counts/interval_lengths
        ratio_data.sort(order='ratio')
        ratio_data = ratio_data[::-1]
        out_stats = self.stats[ratio_data['index']]
        if return_ratio :
            return out_stats, ratio_data
        else :
            return out_stats

    def zero_series_intervals(self) :
        data_size = self.stats.size - 1
        out_data = np.zeros(data_size, self.zero_interval_dtype)
        out_data['end_last'] = self.stats['interval_last_pos'][:-1]
        out_data['start_next'] = self.stats['interval_first_pos'][1:]
        out_data['zero_length'] = out_data['start_next'] - out_data['end_last']
        return out_data
        
        
        
        
        
        
        
        
        
        