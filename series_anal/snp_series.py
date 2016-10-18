# -*- coding: utf-8 -*-
import numpy as np
from ..autosome_snp_data.chrom_snp_series_rdr import chrom_snp_series_factory_cls
from ..autosome_snp_data.chrom_snp_series_data_rdr import chrom_snp_series_data_tables_cls
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
from ..html_display import array_table as html
from itertools import izip


class series_item_class(object) :
    def __init__(self, chrom, series_data, allele_mask) :
        self.chrom = chrom
        self.series_data = series_data
        self.data_index = series_data['data_index']
        self.series_allele_mask = allele_mask
        
    def set_series_snps(self, snp_data, snp_allele_masks) :
        self.snp_data = snp_data
        self.snp_allele_masks = snp_allele_masks

    def struct_data(self) :
        od = self.series_data
        out_data = (od['data_index'], od['first_pos'], od['last_pos'], od['item_count'], od['p90_allele_count'])
        return out_data
        
    def get_snp_pos(self) :
        out_data_size = self.snp_data.size
        out_data = np.zeros(out_data_size, 'u4')
        out_data[:] = self.snp_data['pos']
        return out_data

class snp_series_set(object) :
    total_allele_count = 5008
    cra = country_region_alleles_cls()
    struct_data_dtype = np.dtype([('series_data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'), 
                                 ('snp_count', 'u2'), ('allele_count', 'u2'), ('series_obj', 'O')])
    def __init__(self, chrom, series_data_indexes) :
        self.chrom = chrom
        self.series_data_indexes = series_data_indexes
        self.set_data = None

    def get_series_and_snp_data(self) :
        series_snps_rdr = chrom_snp_series_factory_cls(self.chrom)
        series_data_rdr = chrom_snp_series_data_tables_cls(self.chrom)
        items_series_data, items_allele_masks = series_data_rdr.from_data_indexes(self.series_data_indexes)
        series_data_rdr.close()
        self.series_item_objs = []
        for item_data, item_allele_mask in izip(items_series_data, items_allele_masks) :
            item_series = series_item_class(self.chrom, item_data, item_allele_mask)
            item_snps = series_snps_rdr.item_objs_from_series_data(item_data)
            item_series.set_series_snps(item_snps.snp_data, item_snps.allele_masks)
            self.series_item_objs.append(item_series)    
        series_snps_rdr.close()
        
    def build_struct_data(self) :
        set_data = []
        for item_obj in self.series_item_objs :
            item_data = []
            item_data.extend(item_obj.struct_data())
            item_data.append(item_obj)
            set_data.append(tuple(item_data))
        self.set_data = np.array(set_data, self.struct_data_dtype)

    def region_stats_columns_from_allele_masks(self, allele_masks, num_tag, int_fmt, float_fmt) :
        series_region_counts = self.cra.region_counts_from_allele_masks(allele_masks)
        series_obs_to_pred = self.cra.region_obs_to_pred_from_region_counts(series_region_counts)
        r_c = series_region_counts
        o_p = series_obs_to_pred
        cols = ((('afr',2), ((r_c['afr'], num_tag, int_fmt), (o_p['afr'], num_tag, float_fmt))),
               (('afx',2), ((r_c['afx'], num_tag, int_fmt), (o_p['afx'], num_tag, float_fmt))),
               (('amr',2), ((r_c['amr'], num_tag, int_fmt), (o_p['amr'], num_tag, float_fmt))),
               (('eas',2), ((r_c['eas'], num_tag, int_fmt), (o_p['eas'], num_tag, float_fmt))), 
               (('eur',2), ((r_c['eur'], num_tag, int_fmt), (o_p['eur'], num_tag, float_fmt))),
               (('sas',2), ((r_c['sas'], num_tag, int_fmt), (o_p['sas'], num_tag, float_fmt))),
               (('sax',2), ((r_c['sax'], num_tag, int_fmt), (o_p['sax'], num_tag, float_fmt))))
        return cols                       


    def series_data_html(self, series_data, allele_masks) :
        s_d = series_data
        lengths = s_d['last_pos'] - s_d['first_pos']
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('index', s_d['series_data_index'], num_tag, int_fmt),
                       ('first', s_d['first_pos'], num_tag, big_fmt),
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['snp_count'], num_tag, int_fmt),
                       ('alleles', s_d['allele_count'], num_tag, int_fmt)]
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def snp_data_html(self, snp_data, allele_masks) :
        s_d = snp_data
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('index', s_d['snp_index'], num_tag, int_fmt),
                       ('pos', s_d['pos'], num_tag, big_fmt),
                       ('id', s_d['id']),
                       ('niv', s_d['not_expressed_is_variant']),  
                       ('alleles', s_d['all_count'], num_tag, int_fmt)]
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def series_set_html(self) :
        if self.set_data is None :
            self.get_series_and_snp_data()
            self.build_struct_data()
        out_html = []
        for ind in range(self.set_data.size) :
            bound = ind + 1
            series_item_array = self.set_data[ind:bound]
            series_item = series_item_array[0]
            series_id = str(series_item['snp_count']) + '_' + str(series_item['allele_count'])
            series_item_obj = series_item['series_obj']
            series_item_allele_masks = series_item_obj.series_allele_mask
            series_item_allele_masks.shape = (1, series_item_allele_masks.size)
            out_html.append('<p>')
            out_html.append('<b>' + series_id + ' series</b>')
            out_html.append('</p><p>')
            out_html.append(self.series_data_html(series_item_array, series_item_allele_masks))
            out_html.append('</p>')
            snp_data = series_item_obj.snp_data
            snp_alleles = series_item_obj.snp_allele_masks
            out_html.append('<p>')
            out_html.append('<b>' + series_id + ' series snps</b>')
            out_html.append('</p><p>')
            out_html.append(self.snp_data_html(snp_data, snp_alleles))
        series_html = '\n'.join(out_html)
        return series_html
            
            
            
            
            
            
            
            
            
