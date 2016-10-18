# -*- coding: utf-8 -*-

import numpy as np
from ..autosome_snp_data.chrom_snp_series_data_rdr import chrom_snp_series_data_tables_cls
from ..html_display import array_table as html

class chrom_stats_cls(object) :
    data_dtype = np.dtype([('length', 'u4'), ('snp_count', 'u2')])
    def __init__(self, chrom) :
        self.chrom = chrom
        data_rdr = chrom_snp_series_data_tables_cls(self.chrom)
        series_data = data_rdr.series_data_table[:]
        data_rdr.close()
        series_snp_counts = series_data['item_count']
        m_series = series_snp_counts > 3
        m_not_in_series = np.logical_not(m_series)
        in_series_size = m_series.sum()
        self.in_series_data = np.zeros(in_series_size, self.data_dtype)
        self.in_series_data['snp_count'] = series_snp_counts[m_series]
        self.in_series_data['length'] = series_data['last_pos'][m_series] - series_data['first_pos'][m_series]
        not_in_series_size = m_not_in_series.sum()
        self.not_in_series_data = np.zeros(not_in_series_size, self.data_dtype)
        self.not_in_series_data['snp_count'] = series_snp_counts[m_not_in_series]
        self.not_in_series_data['length'] = series_data['last_pos'][m_not_in_series] - series_data['first_pos'][m_not_in_series]
        
                        
class autosome_snp_series_stats_cls(object) :
    stats_dtype = np.dtype([('chrom', 'S5'), ('series_count', 'u8'), ('series_snp_count', 'u8'), 
                            ('not_in_series_snp_count', 'u8'), ('in_series_ratio', 'f8'), ('series_snp_count_mean', 'u2'), 
                            ('series_snp_count_median', 'u2'), ('series_snp_count_std', 'f8'),
                            ('series_length_mean', 'u4'), ('series_length_median', 'u4'), ('series_length_std', 'f8')])

    classmethod
    def calc_stats(cls, in_series_data, not_in_series_data) :
        series_count = in_series_data.size
        series_snp_count = in_series_data['snp_count'].sum()
        not_in_series_snp_count = not_in_series_data['snp_count'].sum()
        in_series_ratio = float(series_snp_count)/float(series_snp_count+not_in_series_snp_count)
        series_snp_count_mean = in_series_data['snp_count'].mean()
        series_snp_count_median = np.median(in_series_data['snp_count'])
        series_snp_count_std = np.std(in_series_data['snp_count'])
        series_lengths = in_series_data['length']
        series_length_mean = series_lengths.mean()
        series_length_median = np.median(series_lengths)
        series_length_std = np.std(series_lengths)
        out_data = [series_count,  series_snp_count, not_in_series_snp_count, in_series_ratio, 
                    series_snp_count_mean, series_snp_count_median, series_snp_count_std, 
                    series_length_mean, series_length_median, series_length_std]
        return out_data


    def do_stats(self) :
        in_series_chrom_data = []
        not_in_series_chrom_data = []
        self.by_chrom_series_stats = np.zeros(23, self.stats_dtype)
        for chrom in range(1, 23) :
            ind = chrom - 1
            cs_obj = chrom_stats_cls(chrom)
            in_series_data = cs_obj.in_series_data
            not_in_series_data = cs_obj.not_in_series_data
            chrom_stats = [str(chrom)]
            chrom_stats.extend(self.calc_stats(in_series_data, not_in_series_data))
            self.by_chrom_series_stats[ind] = tuple(chrom_stats)
            in_series_chrom_data.append(in_series_data)
            not_in_series_chrom_data.append(not_in_series_data)
        in_series_chrom_data = np.concatenate(in_series_chrom_data)
        not_in_series_chrom_data = np.concatenate(not_in_series_chrom_data)
        chrom_stats = ['all']
        chrom_stats.extend(self.calc_stats(in_series_chrom_data, not_in_series_chrom_data))
        self.by_chrom_series_stats[22] = tuple(chrom_stats)            
            
    def stats_html(self) :
        d = self.by_chrom_series_stats
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('chr', d['chrom']),
                       ('series', d['series_count'], num_tag, big_fmt),
                       ('in_snps', d['series_snp_count'], num_tag, big_fmt),
                       ('not_in_snps', d['not_in_series_snp_count'], num_tag, big_fmt),
                       ('in_ratio', d['in_series_ratio'], num_tag, float_fmt),
                       ('snp_mean', d['series_snp_count_mean'], num_tag, int_fmt),
                       ('snp_med', d['series_snp_count_median'], num_tag, int_fmt),
                       ('snp_std', d['series_snp_count_std'], num_tag, float_fmt),
                       ('len_mean', d['series_length_mean'], num_tag, big_fmt),
                       ('len_med', d['series_length_median'], num_tag, big_fmt),
                       ('len_std', d['series_length_std'].astype('u4'), num_tag, big_fmt)]
        data_html_obj = html.html_table_cls(column_info)                        
        html_table = data_html_obj.assemble_table()
        return html_table
                    
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
                