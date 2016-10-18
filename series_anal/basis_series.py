# -*- coding: utf-8 -*-
import numpy as np
from  ..html_display import array_table as html

class basis_series_cls(object) :
    total_allele_count = 5008
    selection_data_dtype = np.dtype([('data_index', 'u4'), ('allele_count', 'u2'), ('basis_allele_count', 'u2'),
                                     ('series_allele_mask', '?', total_allele_count), ('basis_allele_mask', '?', total_allele_count)])

    def __init__(self, interval_data, analysis_allele_mask=None) :
        self.interval_data = interval_data
        self.analysis_allele_mask = analysis_allele_mask
        self.cra = self.interval_data.cra
        self.interval_first_pos = interval_data.first_pos
        self.interval_last_pos = interval_data.last_pos
        self.series_data = interval_data.series_data
        self.series_allele_masks = interval_data.allele_masks
        if analysis_allele_mask is None :
            self.analysis_allele_mask = np.zeros(self.total_allele_count, '?')
            self.analysis_allele_mask = True
        else :
            test_allele_masks = np.logical_and(self.series_allele_masks, analysis_allele_mask)
            alleles_per_test_mask = test_allele_masks.sum(axis=1)
            m_anal = alleles_per_test_mask > 0
            self.series_data = self.series_data[m_anal]
            self.series_allele_masks = self.series_allele_masks[m_anal]
        self.series_first_pos = self.series_data['first_pos']
        self.series_last_pos = self.series_data['last_pos']
        self.series_lengths = self.series_last_pos - self.series_first_pos + 1
        self.series_snp_counts = self.series_data['item_count']
        self.alleles_per_series = self.series_allele_masks.sum(axis=1)
        self.select_basis_series()

    def select_basis_series(self) :
        selection_data_size = self.series_data.size
        selection_data = np.zeros(selection_data_size, self.selection_data_dtype)
        selection_data['data_index'] = self.series_data['data_index']
        selection_data['series_allele_mask'] = self.series_allele_masks
        selection_data['allele_count'] = self.alleles_per_series
        selection_data.sort(order='allele_count')
        not_selected_allele_mask = self.analysis_allele_mask.copy()
        series_allele_masks = selection_data['series_allele_mask']
        selection_allele_masks = selection_data['basis_allele_mask']
        for index in xrange(selection_data.size) :
            allele_mask = np.logical_and(not_selected_allele_mask, series_allele_masks[index])
            not_selected = np.logical_not(allele_mask)
            not_selected_allele_mask = np.logical_and(not_selected_allele_mask, not_selected)
            selection_allele_masks[index] = allele_mask
        selection_data['basis_allele_count'] = selection_data['basis_allele_mask'].sum(axis=1)
        mgz = selection_data['basis_allele_count'] > 0
        selection_data = selection_data[mgz]
        self.basis_series = selection_data
        self.basis_series.sort(order='data_index')

    def basis_data_from_data_index(self, data_index) :
        ind_data = self.basis_series['data_index'].searchsorted(data_index)
        if ind_data < self.basis_series.size :
            item_data = self.basis_series[ind_data]
            if item_data['data_index'] == data_index :
                return item_data
                
                
    def region_stats_columns_from_allele_masks(self, allele_masks, num_tag, int_fmt, float_fmt, maybe_allele_mask=None) :
        if maybe_allele_mask is None :
            series_region_counts = self.cra.region_counts_from_allele_masks(allele_masks)
        else :
            series_region_counts = self.cra.maybe_region_counts_from_allele_masks(allele_masks, maybe_allele_mask)
        series_obs_to_pred = self.cra.region_obs_to_pred_from_region_counts(series_region_counts, maybe_allele_mask)
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
                
                
    def basis_series_html(self, basis_series_data=None) :
        if basis_series_data is None :
            basis_series_data = self.basis_series
        basis_series_data_indexes = basis_series_data['data_index']
        series_indexes = self.series_data['data_index'].searchsorted(basis_series_data_indexes)
        s_d = self.series_data[series_indexes]
        series_lengths = s_d['last_pos'] - s_d['first_pos']
        basis_allele_counts = basis_series_data['basis_allele_count']
        basis_allele_masks = basis_series_data['basis_allele_mask']
        num_tag = '<td style="text-align: right;">'
        padded_num_tag = '<td style="text-align: right;padding-right: 1em">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('index', s_d['data_index'], num_tag, int_fmt),
                       ('first', s_d['first_pos'], num_tag, big_fmt),
                       ('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', series_lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       ('alleles', s_d['p90_allele_count'], padded_num_tag, int_fmt),
                       ('basis', basis_allele_counts, padded_num_tag, int_fmt)]
        column_info.extend(self.region_stats_columns_from_allele_masks(basis_allele_masks, num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table
                
                
    def html_sorted_by_allele_count(self) :
        sortargs = self.basis_series['allele_count'].argsort()
        sortargs = sortargs[::-1]
        data = self.basis_series[sortargs]
        return self.basis_series_html(data)                        
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                