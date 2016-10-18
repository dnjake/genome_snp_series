# -*- coding: utf-8 -*-
import numpy as np
from autosome_snp_data.chrom_snp_series_data_rdr import chrom_series_in_pos_interval_rdr_cls
from chrom_interval_finder import chrom_interval_zero_pos_cls

class interval_stats_cls(object) :
    total_allele_count = 5008
    grouped_series_dtype = np.dtype([('data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'),
                          ('snp_count', 'u2'), ('allele_count', 'u2'), ('allele_mask', '?', total_allele_count)])
    buckets = np.array([10000, 20000, 40000, 80000, 160000, 320000, 640000], dtype='u4')
    dist_count_dtype = np.dtype([('lower', 'u4'), ('upper', 'u4')])
    stats_dtype = np.dtype([('chrom', 'u2'), ('interval_first_pos', 'u4'), ('interval_last_pos', 'u4'), ('interval_length', 'u4'),
                            ('interval_series_count', 'u2'), ('grouped_snp_count', 'u4'), ('isolated_snp_count', 'u4'),
                            ('max_series_length', 'u4'), ('mean_series_length', 'u4'), ('b0', 'u2'), ('b10', 'u2'),
                            ('b20', 'u2'), ('b40', 'u2'),('b80', 'u2'), ('b160', 'u2'),('b320', 'u2'), ('b640', 'u2'),
                            ('interval_active_allele_count', 'u2'), ('max_active_allele_count', 'u2'), 
                            ('avg_active_allele_count', 'u2')])
                            
    def __init__(self, chrom, interval_first_pos, interval_last_pos) :
        self.chrom = chrom
        self.first_pos = interval_first_pos
        self.last_pos = interval_last_pos
        rdr = chrom_series_in_pos_interval_rdr_cls(self.chrom, self.first_pos, self.last_pos)
        rdr.read_interval_series()
        self.series_data = rdr.series_data
        self.allele_masks = rdr.allele_masks
        self.alleles_per_series = self.allele_masks.sum(axis=1)
        
    def partition_series(self) :
        m = self.series_data['item_count'] > 3
        m_size = m.sum()
        grouped_series_data = np.zeros(m_size, self.grouped_series_dtype)
        for name in ['data_index', 'first_pos', 'last_pos'] :
            grouped_series_data[name] = self.series_data[name][m]
        grouped_series_data['snp_count'] = self.series_data['item_count'][m]
        grouped_series_data['allele_count'] = self.series_data['p90_allele_count'][m]
        grouped_series_data['allele_mask'] = self.allele_masks[m]            
        self.grouped_series_data = grouped_series_data
        nm = np.logical_not(m)
        self.isolated_series_data = self.series_data[nm]
        self.grouped_snp_count = self.grouped_series_data['snp_count'].sum()
        self.isolated_snp_count = self.isolated_series_data['item_count'].sum()
        self.interval_series_count = self.grouped_series_data.size
        
    def calc_series_length_stats(self) :
        series_lengths = self.grouped_series_data['last_pos'] - self.grouped_series_data['first_pos'] + 1        
        self.interval_max_series_length = series_lengths.max()
        self.interval_mean_series_length = series_lengths.mean()
        series_lengths.sort()
        size_dist_indexes = series_lengths.searchsorted(self.buckets)
        dist_counts = np.zeros(self.buckets.size+1, self.dist_count_dtype)
        dist_counts['lower'][1:] = size_dist_indexes
        dist_counts['upper'][:-1] = size_dist_indexes
        dist_counts['upper'][-1] = series_lengths.size
        self.bucket_counts = dist_counts['upper'] - dist_counts['lower']
        
    def alleles_from_pos(self, series_first_pos) :
        m_first = self.grouped_series_data['first_pos'] <= series_first_pos
        m_last = self.grouped_series_data['last_pos'] >= series_first_pos
        m_active = np.logical_and(m_first, m_last)
        active_allele_masks = self.grouped_series_data['allele_mask'][m_active]
        series_per_allele = active_allele_masks.sum(axis=0)
        m_alleles = series_per_allele > 0
        self.sometime_active_alleles[m_alleles] = True
        return m_alleles.sum()

    def calc_allele_stats(self) :
        self.sometime_active_alleles = np.zeros(self.total_allele_count, '?')
        active_allele_count = 0
        interval_pos = self.first_pos
        max_active_allele_count = 0
        length_weighted_active_allele_count = 0
        for item_pos in self.grouped_series_data['first_pos']:
            length_weighted_active_allele_count += active_allele_count*(item_pos - interval_pos)
            interval_pos = item_pos
            active_allele_count = self.alleles_from_pos(item_pos)
            if active_allele_count > max_active_allele_count :
                max_active_allele_count = active_allele_count
        length_weighted_active_allele_count += active_allele_count*(self.last_pos - interval_pos)       
        self.avg_active_allele_count = float(length_weighted_active_allele_count) /float(self.last_pos-self.first_pos)
        self.max_active_allele_count = max_active_allele_count
        self.interval_active_allele_count = self.sometime_active_alleles.sum()
            
    def process_interval(self) :
        self.partition_series()
        self.calc_series_length_stats()
        self.calc_allele_stats()
        
    def stats_for_interval(self) :
        interval_length = self.last_pos - self.first_pos
        out_data = [self.chrom, self.first_pos, self.last_pos, interval_length, self.interval_series_count, 
                    self.grouped_snp_count, self.isolated_snp_count, self.interval_max_series_length,
                    int(self.interval_mean_series_length+0.5)]
        out_data.extend(self.bucket_counts)
        allele_data = [self.interval_active_allele_count, self.max_active_allele_count, int(self.avg_active_allele_count+0.5)]                    
        out_data.extend(allele_data)
        return tuple(out_data)
            
            
class chrom_interval_stats_cls(object) :
    #interval_dtype = np.dtype([('interval_first_pos', 'u4'), ('interval_last_pos', 'u4')])
    def __init__(self, chrom) :
        self.chrom = chrom
        self.intervals = None
        self.stats = None

    def find_intervals(self) :
        zif = chrom_interval_zero_pos_cls(self.chrom)
        self.intervals = zif.series_intervals
        
        '''
        interval_count = zif.zero_snp_series_pos.size
        intervals = np.zeros(interval_count, self.interval_dtype)
        intervals['interval_last_pos'] = zif.zero_snp_series_pos
        intervals['interval_first_pos'][1:] = zif.zero_snp_series_pos[:-1] + 1
        '''        
    def process_intervals(self) :
        if self.intervals is None :
            self.find_intervals()
        out_data = []
        for interval in self.intervals :
            handler = interval_stats_cls(self.chrom, interval['interval_first_pos'], interval['interval_last_pos'])
            #print interval['interval_first_pos'], interval['interval_last_pos'], handler.allele_masks.dtype
            handler.process_interval()
            out_data.append(handler.stats_for_interval())
        self.stats = np.array(out_data, handler.stats_dtype)            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            