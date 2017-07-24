# -*- coding: utf-8 -*-

'''
    snp_data_descr = [ ('snp_index', '<u4'), ('chrom', 'u2'), ('ref', 'S1'), ('alt', 'S1'), ('id', 'S11'), 
                      ('not_expressed_is_variant', 'u1'), ('pos', '<u4'), ('all_count', 'u2') ]


    data_dtype = np.dtype([('data_index', 'u4'), ('first_snp_index', 'u4'), ('chrom', 'u2'), ('first_pos', 'u4'),
                           ('last_pos', 'u4'), ('item_data_start', 'u4'), ('item_count', 'u2'), ('p90_allele_count', 'u2'), 
                           ('afr_obs_to_pred', 'f4'), ('afx_obs_to_pred', 'f4'), ('amr_obs_to_pred', 'f4'), 
                           ('eas_obs_to_pred', 'f4'), ('eur_obs_to_pred', 'f4'), ('sas_obs_to_pred', 'f4'), ('sax_obs_to_pred', 'f4')]  )




'''

import numpy as np
import tables as tb
from genomes_dnj.autosome_snp_data.chrom_snp_series_data_rdr import chrom_series_in_pos_interval_rdr_cls
from genomes_dnj.autosome_snp_data.chrom_snp_data_rdr import chrom_snp_data_tables_cls
from genomes_dnj.autosome_snp_data.chrom_snp_series_rdr import chrom_snp_series_factory_cls

class interval_snp_anal(object) :
    series_item_dtype = np.dtype([('series_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'),
                                  ('snp_count', 'u2'), ('p90_allele_count', 'u2'), ('allele_mask', '?', 5008)])
    snp_item_dtype = np.dtype([('snp_index', 'u4'), ('pos', 'u4'), ('series_index', 'u4'), ('series_snp_count', 'u2'),
                               ('series_allele_count', 'u2'), ('allele_count', 'u2'), ('allele_mask', '?', 5008)])                                  
    chrom = 2        
    first_pos = 135757320L
    last_pos = 136786630L
    
    def __init__(self) :
        rdr = chrom_series_in_pos_interval_rdr_cls(self.chrom, self.first_pos, self.last_pos)
        rdr.read_interval_series()
        series_data = rdr.series_data
        print 'series_data_size', series_data.size
        self.snp_series = np.zeros(series_data.size, self.series_item_dtype)
        in_names = ('data_index', 'first_pos', 'last_pos', 'item_count', 'p90_allele_count')
        out_names = ('series_index', 'first_pos', 'last_pos', 'snp_count', 'p90_allele_count')
        data_names = zip(in_names, out_names)
        for in_name, out_name in data_names:
            self.snp_series[out_name] = series_data[in_name]
        self.snp_series['allele_mask'] = rdr.allele_masks
        self.read_snps()
        self.add_series_data(series_data)
        
    def read_snps(self) :
        st = chrom_snp_data_tables_cls(self.chrom)
        range_start_pos = self.first_pos
        range_bound_pos = self.last_pos       
        gcond = '(pos >= ' + str(range_start_pos) + ')'
        lcond = '(pos <= ' + str(range_bound_pos) + ')'
        range_condition = gcond + ' & ' + lcond
        in_snp_data = st.snp_data_table.read_where(range_condition)
        print 'in_snp_data', in_snp_data.size
        self.snp_data = np.zeros(in_snp_data.size, self.snp_item_dtype)
        in_names = ('snp_index', 'pos', 'all_count')
        out_names = ('snp_index', 'pos', 'allele_count')
        names = zip(in_names, out_names)
        for in_name, out_name in names :
            self.snp_data[out_name] = in_snp_data[in_name]
        data_indexes = self.snp_data['snp_index']
        bp_alleles = st.snp_bitpacked_allele_values_table[data_indexes]        
        self.snp_data['allele_mask'] = np.unpackbits(bp_alleles['bitpacked_values'], axis=1)
        st.close()
        
    def add_series_data(self, series_data) :
        ss = chrom_snp_series_factory_cls(self.chrom)
        for series_item in series_data :
            series_index = series_item['data_index']
            snp_count = series_item['item_count']
            series_allele_count = series_item['p90_allele_count']
            series_snp_indexes = ss.item_snp_indexes_from_series_data(series_item)
            series_snp_indexes.sort()
            assert snp_count == series_snp_indexes.size
            snp_data_indexes = self.snp_data['snp_index'].searchsorted(series_snp_indexes)
            self.snp_data['series_index'][snp_data_indexes] = series_index
            self.snp_data['series_snp_count'][snp_data_indexes] = snp_count
            self.snp_data['series_allele_count'][snp_data_indexes] = series_allele_count
        ss.close()
            
        
data = interval_snp_anal()        

filters = tb.Filters(complib='zlib', complevel=5)
h5_out = tb.open_file('lct_interval_snps.h5', 'w', filers=filters)
out_series = h5_out.create_table('/', 'snp_series_data', description=data.snp_series.dtype)
out_series.append(data.snp_series)
out_snps = h5_out.create_table('/', 'snp_data', description=data.snp_data.dtype)
out_snps.append(data.snp_data)
h5_out.close()
print 'done'


        
        
        
        
        
        
        
        
        
        
        
        