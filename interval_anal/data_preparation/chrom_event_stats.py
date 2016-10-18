# -*- coding: utf-8 -*-

import numpy as np
from autosome_snp_data.chrom_snp_series_data_rdr import chrom_series_in_pos_interval_rdr_cls
from chrom_interval_finder import chrom_interval_zero_pos_cls

class interval_event_stats_cls(object) :
    total_allele_count = 5008
    event_dtype = np.dtype([('event_pos', 'u4'), ('is_first_pos', '?')])
    event_stats_dtype = np.dtype([('event_pos', 'u4'), ('is_first_pos', '?'), ('active_series_avg_length', 'u4'),
                                  ('active_series', 'u2'), ('active_series_g_100', 'u2'), ('active_series_g_500', 'u2'), 
                                  ('active_series_snp_count', 'u2'), ('active_series_allele_count', 'u2')])
    grouped_series_dtype = np.dtype([('data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'),
                          ('snp_count', 'u2'), ('allele_count', 'u2'), ('allele_mask', '?', total_allele_count)])
                            
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
        event_size = grouped_series_data.size
        grouped_series_events = np.zeros(2*event_size, self.event_dtype)
        grouped_series_events['event_pos'][:event_size] = grouped_series_data['first_pos']
        grouped_series_events['is_first_pos'][:event_size] = True
        grouped_series_events['event_pos'][event_size:] = grouped_series_data['last_pos']
        grouped_series_events.sort(order='event_pos')
        self.grouped_series_events = grouped_series_events
        
        
    def active_stats_from_pos(self, pos) :
        m_first = self.grouped_series_data['first_pos'] <= pos
        m_last = self.grouped_series_data['last_pos'] > pos
        m_active = np.logical_and(m_first, m_last)
        active_series = self.grouped_series_data[m_active]
        active_series_count = active_series.size
        active_series_avg_length = 0
        active_series_g_100 = 0
        active_series_g_500 = 0
        if active_series_count > 0 :
            active_series_lengths = active_series['last_pos'] - active_series['first_pos']
            active_series_avg_length = active_series_lengths.mean()
            m_100 = active_series_lengths >= 100000
            active_series_g_100 = m_100.sum()
            m_500 = active_series_lengths >= 500000
            active_series_g_500 = m_500.sum()
        active_snp_count = active_series['snp_count'].sum()
        active_allele_masks = self.grouped_series_data['allele_mask'][m_active]
        series_per_allele = active_allele_masks.sum(axis=0)
        m_alleles = series_per_allele > 0
        active_allele_count = m_alleles.sum()
        out_data = (active_series_avg_length, active_series_count, active_series_g_100, 
                    active_series_g_500, active_snp_count, active_allele_count)
        return out_data
                    
    def do_event_stats(self) :
        out_data = []
        for event_pos, is_first_pos in self.grouped_series_events:
            event_data = [event_pos, is_first_pos]
            event_data.extend(self.active_stats_from_pos(event_pos))
            out_data.append(tuple(event_data))
        out_data = np.array(out_data, self.event_stats_dtype)
        return out_data
            
    def process_interval(self) :
        self.partition_series()
        interval_data = self.do_event_stats()
        return interval_data



class chrom_event_stats_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.intervals = None
        self.stats = None

    def find_intervals(self) :
        zif = chrom_interval_zero_pos_cls(self.chrom)
        self.intervals = zif.series_intervals
        
    def process_intervals(self) :
        if self.intervals is None :
            self.find_intervals()
        out_data = []
        for interval in self.intervals :
            handler = interval_event_stats_cls(self.chrom, interval['interval_first_pos'], interval['interval_last_pos'])
            #print interval['interval_first_pos'], interval['interval_last_pos'], handler.allele_masks.dtype
            out_data.append(handler.process_interval())
        self.stats = np.concatenate(out_data)            

















