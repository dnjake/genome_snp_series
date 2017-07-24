# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
import os
np.set_printoptions(threshold=5000)
from genomes_dnj.lct_interval.lct_interval_plot_context import plot_context_cls

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

'''
>>> snp_data.dtype
dtype([('snp_index', '<u4'), ('pos', '<u4'), ('series_index', '<u4'), ('series_snp_count', '<u2'), 
('series_allele_count', '<u2'), ('allele_count', '<u2'), ('allele_mask', '?', (5008,))])

>>> self.series_data.dtype
dtype([('series_index', '<u4'), ('first_pos', '<u4'), ('last_pos', '<u4'), ('snp_count', '<u2'), 
('p90_allele_count', '<u2'), ('allele_mask', '?', (5008,))])

'''

class snp_and_series_cls(object) :
    file_name = 'lct_interval_snps.h5'
    def __init__(self) :
        file_path = os.path.join(mod_dir, self.file_name)
        h5 = tb.open_file(file_path, 'r')
        self.snp_data = h5.root.snp_data[:]
        self.series_data = h5.root.snp_series_data[:]
        h5.close()
        self.snp_indexes = self.snp_data['snp_index']
        self.series_indexes = self.series_data['series_index']

    def series_from_index(self, series_index) :
        data_index = self.series_indexes.searchsorted(series_index)
        if series_index == self.series_indexes[data_index] :
            return self.series_data[data_index]
                
    def series_allele_mask_from_index(self, series_index) :
        series_item = self.series_from_index(series_index)
        return series_item['allele_mask']

    def series_snps_from_index(self, series_index) :
        m = self.snp_data['series_index'] == series_index
        return self.snp_data[m]

'''
what I want to know is series snps expressed by allees that
express the series and series snps by alleles that don't express the series
the real question is whether alleles that express the same count of snps
are expressing the same snps

another question is whether in series alleles are expressing ungrouped
snps within the series extent or snps from other series within that extent
'''
def pos_pair_from_snp_indexes(index_data) :
    out_data = []
    for obj, obj_label, low_index in index_data:
        high_index = low_index + 1
        obj_pos = obj.series_pos
        low_pos = obj_pos[low_index]
        high_pos = obj_pos[high_index]
        out_pos = (obj_label, low_pos, high_pos)
        out_data.append(out_pos)
    return out_data

pos_data_dtype = np.dtype([('label', 'S20'), ('low', 'u4'), ('high', 'u4')])    
def array_from_pos_pair_data(data) : 
    array_size = len(data) + 1
    data_array = np.zeros(array_size, pos_data_dtype)
    work_data = data_array[:-1]
    work_data[:] = data
    work_low = work_data['low']
    low_max = work_low.max()
    work_high = work_data['high']
    high_min = work_high.min()
    data_last = data_array[-1]
    data_last['label'] = 'min max'
    data_last['low'] = low_max
    data_last['high'] = high_min
    return data_array

def pair_array_from_indexes(index_data) :
    pair_data = pos_pair_from_snp_indexes(index_data)
    array_data = array_from_pos_pair_data(pair_data)
    return array_data

    
class series_anal_cls(object) :
    unique_count_dtype = np.dtype([('count', 'u2'), ('snps', 'u2')])
       
    def __init__(self, series_index, data=None) :
        self.series_index = series_index        
        self.data = data
        if self.data is None :
            self.data = snp_and_series_cls()
        self.series = self.data.series_from_index(self.series_index)
        self.series_snps = self.data.series_snps_from_index(self.series_index)
        self.series_pos = self.series_snps['pos']

    def unique_snps_per_allele(self, allele_mask=None) :
        if allele_mask is None :
            allele_mask = np.ones(5008, '?')
        snp_allele_masks = np.logical_and(self.series_snps['allele_mask'], allele_mask)
        snps_per_allele = snp_allele_masks.sum(axis=0)
        snps_per_allele = snps_per_allele[allele_mask]
        snp_counts, instance_counts = np.unique(snps_per_allele, return_counts=True)
        data = zip(snp_counts, instance_counts)
        data_a = np.array(data, dtype=self.unique_count_dtype)
        return data_a
                
    def snps_from_aps_value(self, value, allele_mask=None) :
        snp_allele_masks = self.series_snps['allele_mask']
        aps = snp_allele_masks.sum(axis=0)
        m = aps == value
        if allele_mask is not None :  
            m = np.logical_and(m, allele_mask)
        snp_allele_masks = np.logical_and(snp_allele_masks, m)
        spa = snp_allele_masks.sum(axis=1)
        return spa, m
            
    def snp_allele_mask_from_index(self, index, allele_mask=None) :
        snp_allele_masks = self.series_snps['allele_mask']
        if allele_mask is not None :  
            snp_allele_masks = np.logical_and(snp_allele_masks, allele_mask)
        return snp_allele_masks[index]

    def snp_allele_counts_and_indexes(self, allele_mask) :
        snp_allele_masks = np.logical_and(self.series_snps['allele_mask'], allele_mask)
        snps_per_allele = snp_allele_masks.sum(axis=0)
        m = snps_per_allele > 0
        snps_per_allele = snps_per_allele[m]
        allele_indexes = np.where(m)[0]
        return snps_per_allele, allele_indexes

    def alleles_per_snp(self, allele_mask) :
        snp_allele_masks = np.logical_and(self.series_snps['allele_mask'], allele_mask)
        alleles_per_snp = snp_allele_masks.sum(axis=1)
        return alleles_per_snp
        
    def get_in_series_snps_per_allele(self) :
        allele_mask = self.series['allele_mask']
        counts_indexes = self.snp_allele_counts_and_indexs(allele_mask)
        self.in_series_snps_per_allele, self.in_series_allele_indexes = counts_indexes
        
    def get_out_series_snps_per_allele(self) :
        allele_mask = self.series['allele_mask']
        allele_mask = np.logical_not(allele_mask)
        counts_indexes = self.snp_allele_counts_and_indexs(allele_mask)
        self.out_series_snps_per_allele, self.out_series_allele_indexes = counts_indexes

    def snps_per_allele(self, snp_mask) :
        snp_allele_masks = self.series_snps['allele_mask']
        snps_per_allele = snp_allele_masks[snp_mask].sum(axis=0)
        return snps_per_allele
        
'''
Need to select snps in fragment and alleles that express it
find the allele mask in a static
'''

class series_fragment_cls(object) :
    def __init__(self, series, fragment_snp_mask,  fragment_allele_mask) :
        self.series = series
        self.fragment_snp_mask = fragment_snp_mask
        self.fragment_allele_mask = fragment_allele_mask


'''
find snps_allele
find max count
find snps expressed by min_match of alleles that express the max count
find alleles that express min match of those snps
'''
class series_fragment_finder_cls(object) :        
    aps_dtype = np.dtype([('index', 'u2'), ('snp_count', 'u2')])
    with_pos_dtype = np.dtype([('index', 'u2'), ('pos', 'u4'), ('snp_count', 'u2')])
    def __init__(self, series, allele_mask) :
        self.series = series
        self.allele_mask = allele_mask
        self.count_data = series.unique_snps_per_allele(self.allele_mask)
        self.arg_max_snps = self.count_data['snps'].argmax()
        self.max_count = self.count_data['count'][self.arg_max_snps]
        max_snps_alleles_per_snp, self.max_snps_allele_mask = self.series.snps_from_aps_value(self.max_count, self.allele_mask)
        self.max_snps_alleles_per_snp = np.zeros(max_snps_alleles_per_snp.size, self.aps_dtype)
        self.max_snps_alleles_per_snp['snp_count'] = max_snps_alleles_per_snp
        self.max_snps_alleles_per_snp['index'] = np.arange(self.max_snps_alleles_per_snp.size, dtype='u2')
    
    def fragment_from_max_snps(self, min_match=0.9) :
        allele_count = self.max_snps_allele_mask.sum()
        fragment_snp_mask = self.max_snps_alleles_per_snp >= min_match*allele_count
        snps_per_allele = self.series.snps_per_allele(fragment_snp_mask)
        allele_mask = snps_per_allele >= min_match*self.max_count
        allele_mask = np.logical_and(allele_mask, self.allele_mask)
        fragment = series_fragment_cls(self.series, fragment_snp_mask, allele_mask)
        return fragment
        
    def alleles_per_snp_from_count(self, count) :
        alleles_per_snp, allele_mask = self.series.snps_from_aps_value(count, self.allele_mask)
        data = np.zeros(alleles_per_snp.size, self.aps_dtype)
        data['snp_count'] = alleles_per_snp
        data['index'] = np.arange(data.size, dtype='u2')
        return data
                
    def aps_with_pos(self, count_data) :
        d = np.zeros(count_data.size, self.with_pos_dtype)
        for name in ['index', 'snp_count'] :
            d[name] = count_data[name]
        d['pos'] = self.series.series_snps['pos']
        return d
        
    def am_from_index(self, index, count_allele_mask=None) :
        am = self.series.series_snps['allele_mask'][index]
        if count_allele_mask is None :
            count_allele_mask = self.allele_mask
        am = np.logical_and(am, count_allele_mask)
        return am

    def am_from_spa_count_and_index(self, spa_count, index) :
        alleles_per_snp, allele_mask = self.series.snps_from_aps_value(spa_count, self.allele_mask)
        index_allele_mask = self.series.series_snps['allele_mask'][index]
        allele_mask = np.logical_and(allele_mask, index_allele_mask)
        return allele_mask




       