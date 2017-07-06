# -*- coding: utf-8 -*-

import numpy as np
from ..autosome_snp_data.chrom_snp_series_data_rdr import chrom_series_in_pos_interval_rdr_cls
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
from ..html_display import array_table as html
#cra = country_region_alleles_cls()

class interval_data_cls(object) :
    cra = country_region_alleles_cls()
    total_allele_count = 5008
    min_series_snp_count = 4
    to_layout_dtype = np.dtype([('data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'), ('allele_count', 'u2')])
    @classmethod
    def to_layout_series_data(cls, p90_series_data) :
        out_data_dtype = cls.to_layout_dtype
        layout_series_data = np.zeros(p90_series_data.size, out_data_dtype)
        for name in ['data_index', 'first_pos', 'last_pos'] :
            layout_series_data[name] = p90_series_data[name]
        layout_series_data['allele_count'] = p90_series_data['p90_allele_count']
        return layout_series_data
        
    def __init__(self, chrom, interval_first_pos, interval_last_pos) :
        self.chrom = chrom
        self.first_pos = interval_first_pos
        self.last_pos = interval_last_pos
        rdr = chrom_series_in_pos_interval_rdr_cls(self.chrom, self.first_pos, self.last_pos)
        rdr.read_interval_series()
        m = rdr.series_data['item_count'] >= self.min_series_snp_count
        self.series_data = rdr.series_data[m]
        self.data_indexes = self.series_data['data_index']
        self.allele_masks = rdr.allele_masks[m]
        self.alleles_per_series = self.allele_masks.sum(axis=1)

    def or_allele_mask(self, data_indexes) :
        interval_indexes = self.series_data['data_index'].searchsorted(data_indexes)
        out_mask = np.zeros(self.total_allele_count, '?')
        for index in interval_indexes :
            index_mask = self.allele_masks[index]
            out_mask = np.logical_or(out_mask, index_mask)
        return out_mask
        
    def and_allele_mask(self, data_indexes) :
        interval_indexes = self.series_data['data_index'].searchsorted(data_indexes)
        out_mask = self.allele_masks[interval_indexes[0]].copy()
        for index in interval_indexes[1:] :                    
            index_mask = self.allele_masks[index]
            out_mask = np.logical_and(out_mask, index_mask)
        return out_mask

    def or_and_not_allele_mask(self, or_data_indexes, not_or_data_indexes) :
        or_allele_mask = self.or_allele_mask(or_data_indexes)
        not_or_allele_mask = self.or_allele_mask(not_or_data_indexes)
        not_or_allele_mask = np.logical_not(not_or_allele_mask)
        return np.logical_and(or_allele_mask, not_or_allele_mask)

    def and_and_not_or_allele_mask(self, and_data_indexes, not_or_data_indexes) :
        and_allele_mask = self.and_allele_mask(and_data_indexes)
        not_or_allele_mask = self.or_allele_mask(not_or_data_indexes)
        not_or_allele_mask = np.logical_not(not_or_allele_mask)
        return np.logical_and(and_allele_mask, not_or_allele_mask)

    def not_or_allele_mask(self, not_data_indexes) :
        allele_mask = self.or_allele_mask(not_data_indexes)
        return np.logical_not(allele_mask)
        
    def not_and_allele_mask(self, not_data_indexes) :
        allele_mask = self.and_allele_mask(not_data_indexes)
        return np.logical_not(allele_mask)
        
    def yes_no_data_indexes(self, yes_data_indexes=None, no_data_indexes=None) :
        if yes_data_indexes is None :
            result = self.not_or_allele_mask(no_data_indexes)
        elif no_data_indexes is None :
            result = self.and_allele_mask(yes_data_indexes)
        else :
            result = self.and_and_not_or_allele_mask(yes_data_indexes, no_data_indexes)
        result = result.astype('?')
        return result

    def yes_allele_mask_or_indexes(self, yes_allele_mask, data_indexes) :
        or_allele_mask = self.or_allele_mask(data_indexes)
        return np.logical_and(yes_allele_mask, or_allele_mask)
        
    def yes_allele_mask_and_indexes(self, yes_allele_mask, data_indexes) :
        and_allele_mask = self.and_allele_mask(data_indexes)
        return np.logical_and(yes_allele_mask, and_allele_mask)
 
    def yes_allele_mask_not_or_indexes(self, yes_allele_mask, no_data_indexes) :
        no_allele_mask = self.not_or_allele_mask(no_data_indexes)
        return np.logical_and(yes_allele_mask, no_allele_mask)
        
    def yes_allele_mask_not_and_indexes(self, yes_allele_mask, no_data_indexes) :
        no_allele_mask = self.not_and_allele_mask(no_data_indexes)
        return np.logical_and(yes_allele_mask, no_allele_mask)

    def allele_mask_from_data_index(self, data_index) :
        interval_index = self.series_data['data_index'].searchsorted(data_index)
        out_mask = self.allele_masks[interval_index].copy()
        out_mask = out_mask.astype('?')
        return out_mask
        
    def superset_series_mask_from_allele_indexes(self, allele_indexes, min_match=0.9) :
        test_allele_masks = self.allele_masks[:,allele_indexes]
        match_alleles_per_series = test_allele_masks.sum(axis=1)
        match_thresh = min_match*float(allele_indexes.size)
        mask = match_alleles_per_series >= match_thresh
        return mask, match_alleles_per_series
        
    def subset_series_mask_from_allele_indexes(self, allele_indexes, min_match=0.9) :
        test_allele_masks = self.allele_masks[:,allele_indexes]
        match_alleles_per_series = test_allele_masks.sum(axis=1)
        mask = match_alleles_per_series >= min_match*self.alleles_per_series.astype('f4')
        return mask, match_alleles_per_series

    def superset_data_from_data_index(self, data_index, min_match=0.9) :
        index = self.data_indexes.searchsorted(data_index)
        index_allele_mask = self.allele_masks[index]
        index_allele_indexes = np.where(index_allele_mask)[0]
        data_mask, match_alleles_per_series = self.superset_series_mask_from_allele_indexes(index_allele_indexes, min_match)
        #data_mask[index] = False
        series_data = self.series_data[data_mask]
        allele_data = self.allele_masks[data_mask]
        match_alleles_per_series = match_alleles_per_series[data_mask]
        return series_data, allele_data, match_alleles_per_series
        
    def subset_data_from_allele_mask(self, allele_mask) :        
        allele_indexes = np.where(allele_mask)[0]
        data_mask, match_alleles_per_series = self.subset_series_mask_from_allele_indexes(allele_indexes)
        series_data = self.series_data[data_mask]
        allele_data = self.allele_masks[data_mask]
        match_alleles_per_series = match_alleles_per_series[data_mask]
        return series_data, allele_data, match_alleles_per_series

    def subset_data_from_data_index(self, data_index) :
        index = self.data_indexes.searchsorted(data_index)
        index_allele_mask = self.allele_masks[index]
        return self.subset_data_from_allele_mask(index_allele_mask)

                
    def superset_data_from_match_allele_mask(self, match_allele_mask, min_match=0.9) :
        #match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        match_allele_count = match_allele_mask.sum()
        if match_allele_count == 0 :
            matched_series = np.empty(0, self.series_data.dtype)
            matched_series_allele_masks = np.empty((0, self.total_allele_count), '?')
            series_match_masks = np.empty((0, self.total_allele_count), '?')
        else :
            series_match_masks  = np.logical_and(self.allele_masks, match_allele_mask)
            match_alleles_per_series = series_match_masks.sum(axis=1)
            min_match_count = min_match*float(match_allele_count)
            match_mask = match_alleles_per_series >= min_match_count
            matched_series = self.series_data[match_mask]
            matched_series_allele_masks = self.allele_masks[match_mask]
            series_match_masks = series_match_masks[match_mask]
        return matched_series, matched_series_allele_masks, series_match_masks

    def superset_data_from_yes_no_indexes(self, yes_indexes, no_indexes=None, min_match=0.9) :
        match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        matched_series, matched_series_allele_masks, series_match_masks = (
                                      self.superset_data_from_match_allele_mask(match_allele_mask, min_match))
        return match_allele_mask, matched_series, matched_series_allele_masks, series_match_masks
        

    def subset_data_from_match_allele_mask(self, match_allele_mask, min_match=0.9) :
        series_match_masks  = np.logical_and(self.allele_masks, match_allele_mask)
        match_alleles_per_series = series_match_masks.sum(axis=1)
        min_match_counts = min_match*self.alleles_per_series.astype('f4')
        match_mask = match_alleles_per_series >= min_match_counts
        matched_series = self.series_data[match_mask]
        matched_series_allele_masks = self.allele_masks[match_mask]
        series_match_masks = series_match_masks[match_mask]
        return matched_series, matched_series_allele_masks, series_match_masks


    def subset_data_from_yes_no_indexes(self, yes_indexes, no_indexes=None, min_match=0.9) :
        match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        return self.subset_data_from_match_allele_mask(match_allele_mask, min_match)
        

    def series_data_from_data_indexes(self, data_indexes) :
        indexes = self.data_indexes.searchsorted(data_indexes)
        series_data = self.series_data[indexes]
        allele_masks = self.allele_masks[indexes]
        return series_data, allele_masks

    def expressed_series_from_allele_mask(self, allele_mask, min_match=0.0001) :
        series_match_masks  = np.logical_and(self.allele_masks, allele_mask)
        matches_per_series = series_match_masks.sum(axis=1)
        match_ratios = matches_per_series.astype('f4')/self.alleles_per_series.astype('f4')
        m = match_ratios > min_match
        matched_series = self.series_data[m]
        matched_series_allele_masks = self.allele_masks[m]
        series_match_masks = series_match_masks[m]
        return matched_series, matched_series_allele_masks, series_match_masks
        
    def series_with_allele_mask_matches(self, allele_mask, min_match=0.9) :
        series_match_masks  = np.logical_and(self.allele_masks, allele_mask)
        matches_per_series = series_match_masks.sum(axis=1)
        match_ratios = matches_per_series.astype('f4')/float(allele_mask.sum())
        m = match_ratios > min_match
        matched_series = self.series_data[m]
        matched_series_allele_masks = self.allele_masks[m]
        series_match_masks = series_match_masks[m]
        return matched_series, matched_series_allele_masks, series_match_masks

    
    def maybe_series_data(self, data_indexes, maybe_allele_mask):
        series_data, allele_masks = self.series_data_from_data_indexes(data_indexes)
        allele_masks = np.logical_and(allele_masks, maybe_allele_mask)
        return series_data, allele_masks

    def maybe_allele_mask_from_no_data_indexes(self, no_data_indexes) :
        or_allele_mask = self.or_allele_mask(no_data_indexes) 
        maybe_allele_mask = np.logical_not(or_allele_mask)
        return maybe_allele_mask

    def common_series_statistics(self, yes_data_indexes, maybe_allele_mask) :
        series_data, allele_masks = self.series_data_from_data_indexes(yes_data_indexes)
        allele_masks = np.logical_and(allele_masks, maybe_allele_mask)
        series_count = series_data.size
        series_per_allele = allele_masks.sum(axis=0)
        common_allele_mask = series_per_allele == series_count
        common_snp_count = series_data['item_count'].sum()
        common_first_pos = series_data['first_pos'].min()
        common_last_pos = series_data['last_pos'].max()
        common_series_length = common_last_pos - common_first_pos
        common_allele_count = common_allele_mask.sum()
        out_data = (common_first_pos, common_last_pos, common_series_length,
                    common_snp_count, common_allele_count, common_allele_mask)
        return out_data
                    
    def region_counts_from_allele_masks(self, allele_masks) :
        return self.cra.region_counts_from_allele_masks(allele_masks)

    def series_html(self, series_data, num_tag, int_fmt, big_fmt) :
        s_d = series_data
        lengths = s_d['last_pos'] - s_d['first_pos']
        #column_info = (('index', s_d['data_index'], num_tag, int_fmt),
        column_info = [('first', s_d['first_pos'], num_tag, big_fmt),
                       ('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       ('alleles', s_d['p90_allele_count'], num_tag, int_fmt)]
        return column_info

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
        
    def maybe_series_html(self, data_indexes, maybe_allele_mask):
        series_data, allele_masks = self.maybe_series_data(data_indexes, maybe_allele_mask)
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = self.series_html(series_data, num_tag, int_fmt, big_fmt)
        allele_masks = np.logical_and(allele_masks, maybe_allele_mask)
        maybe_match_counts = allele_masks.sum(axis=1)
        column_info.append(('matches', maybe_match_counts, num_tag, int_fmt))
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt, maybe_allele_mask))
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def series_data_html(self, data_indexes) :
        series_data, allele_masks = self.series_data_from_data_indexes(data_indexes)
        s_d = series_data
        lengths = s_d['last_pos'] - s_d['first_pos']
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        #column_info = (('index', s_d['data_index'], num_tag, int_fmt),
        column_info = [('first', s_d['first_pos'], num_tag, big_fmt),
                       ('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       ('alleles', s_d['p90_allele_count'], num_tag, int_fmt)]
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table
        
    def match_set_html(self, series_data, allele_masks, match_allele_counts) :
        sort_indexes = np.argsort(series_data, order='p90_allele_count')
        sort_indexes = sort_indexes[::-1]
        series_data = series_data[sort_indexes]
        allele_masks = allele_masks[sort_indexes]
        match_allele_counts = match_allele_counts[sort_indexes]
        s_d = series_data
        allele_counts = s_d['p90_allele_count']
        allele_counts_f = allele_counts.astype('f4')
        match_allele_counts_f = match_allele_counts.astype('f4')
        match_ratios = np.zeros(allele_counts.size, dtype='f4')
        m = allele_counts_f > 0.0
        match_ratios[m] = match_allele_counts_f[m]/allele_counts_f[m]
        lengths = s_d['last_pos'] - s_d['first_pos']
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        #column_info = (('index', s_d['data_index'], num_tag, int_fmt),
        column_info = [('first', s_d['first_pos'], num_tag, big_fmt),
                       ('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       ('alleles', s_d['p90_allele_count'], num_tag, int_fmt),
                       (('matches',2), ((match_allele_counts, num_tag, int_fmt), (match_ratios, num_tag, float_fmt)))]
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def superset_html_from_data_index(self, data_index, min_match=0.9) :
        series_data, allele_masks, match_allele_counts = self.superset_data_from_data_index(data_index, min_match)
        return self.match_set_html(series_data, allele_masks, match_allele_counts)

    '''
    The series allele masks has the alleles that express the series
    The match allele masks have the alleles that match the particular test
    '''        
    '''
    def superset_html(self, series_data, series_allele_masks, match_allele_masks) :                               
        sort_indexes = np.argsort(series_data, order='p90_allele_count')
        sort_indexes = sort_indexes[::-1]
        series_data = series_data[sort_indexes]
        series_allele_masks = series_allele_masks[sort_indexes]
        match_allele_masks = match_allele_masks[sort_indexes]
        s_d = series_data
        lengths = s_d['last_pos'] - s_d['first_pos']
        match_alleles_per_series = match_allele_masks.sum(axis=1)
        match_alles_per_series_f = match_alleles_per_series.astype('f4')
        alleles_per_series_f = s_d['p90_allele_count'].astype('f4')
        allele_ratios = match_alles_per_series_f/alleles_per_series_f
        #match_allele_counts_f = match_allele_counts.astype('f4')
        #match_ratios = np.zeros(allele_counts.size, dtype='f4')
        #m = allele_counts_f > 0.0
        #match_ratios[m] = match_allele_counts_f[m]/allele_counts_f[m]
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('index', s_d['data_index'], num_tag, int_fmt),
                       ('first', s_d['first_pos'], num_tag, big_fmt),
                       #('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       (('alleles',2), ((s_d['p90_allele_count'], num_tag, int_fmt), (allele_ratios, num_tag, float_fmt))),
                       (('matches',2), ((match_alleles_per_series, num_tag, int_fmt), (allele_ratios, num_tag, float_fmt)))]
        column_info.extend(self.region_stats_columns_from_allele_masks(match_allele_masks, 
                                                               num_tag, int_fmt, float_fmt))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table
        '''
        
    def superset_html_from_yes_no_indexes(self, yes_indexes, no_indexes=None, min_match=0.9) :
        data = self.superset_data_from_yes_no_indexes(yes_indexes, no_indexes, min_match)
        focus_allele_mask, matched_series, matched_series_allele_masks, series_match_masks = data
        return self.match_mask_html(focus_allele_mask, matched_series, matched_series_allele_masks, series_match_masks)

    def superset_html_from_mask(self, allele_mask, min_match=0.9) :
        data = self.superset_data_from_match_allele_mask(allele_mask, min_match)
        matched_series, matched_series_allele_masks, series_match_masks = data
        return self.match_mask_html(allele_mask, matched_series, matched_series_allele_masks, series_match_masks)

    def subset_html_from_data_index(self, data_index) :
        series_data, allele_masks, match_allele_counts = self.subset_data_from_data_index(data_index)
        return self.match_set_html(series_data, allele_masks, match_allele_counts)
        

    '''
    problem is that I am not accounting for the removal of the no allele indexes
    '''
    def common_series_html(self, yes_data_indexes, maybe_allele_mask) :        
        stats = self.common_series_statistics(yes_data_indexes, maybe_allele_mask)
        first_pos, last_pos, length, snp_count, allele_count, allele_mask = stats
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('first', [first_pos], num_tag, big_fmt),
                       ('last', [last_pos], num_tag, big_fmt), 
                       ('length', [length], num_tag, big_fmt),
                       ('snps', [snp_count], num_tag, int_fmt),
                       ('alleles', [allele_count], num_tag, int_fmt)]
        #print allele_mask.sum()
        allele_masks = allele_mask.reshape([1, self.total_allele_count])
        column_info.extend(self.region_stats_columns_from_allele_masks(allele_masks, num_tag, int_fmt, float_fmt, maybe_allele_mask))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table



    def match_mask_html(self, focus_allele_mask, series_data, series_allele_masks, match_allele_masks, maybe_allele_mask=None) :                               
        sort_indexes = np.argsort(series_data, order='p90_allele_count')
        sort_indexes = sort_indexes[::-1]
        series_data = series_data[sort_indexes]
        series_allele_masks = series_allele_masks[sort_indexes]
        match_allele_masks = match_allele_masks[sort_indexes]
        s_d = series_data
        lengths = s_d['last_pos'] - s_d['first_pos']
        match_alleles_per_series = match_allele_masks.sum(axis=1)
        total_focus_alleles = float(focus_allele_mask.sum())
        match_alles_per_series_f = match_alleles_per_series.astype('f4')
        match_ratios = match_alles_per_series_f/total_focus_alleles
        alleles_per_series_f = s_d['p90_allele_count'].astype('f4')
        allele_ratios = match_alles_per_series_f/alleles_per_series_f
        num_tag = '<td style="text-align: right;">'
        int_fmt =  '{:d}'
        big_fmt =  '{:,}'
        float_fmt = '{:.2f}'
        column_info = [('index', s_d['data_index'], num_tag, int_fmt),
                       ('first', s_d['first_pos'], num_tag, big_fmt),
                       #('last', s_d['last_pos'], num_tag, big_fmt), 
                       ('length', lengths, num_tag, big_fmt),
                       ('snps', s_d['item_count'], num_tag, int_fmt),
                       (('alleles',2), ((s_d['p90_allele_count'], num_tag, int_fmt), (allele_ratios, num_tag, float_fmt))),
                       (('matches',2), ((match_alleles_per_series, num_tag, int_fmt), (match_ratios, num_tag, float_fmt)))]
        column_info.extend(self.region_stats_columns_from_allele_masks(match_allele_masks, 
                                                               num_tag, int_fmt, float_fmt, maybe_allele_mask))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def match_series_html(self, allele_mask, min_match=0.0001) :
        matched_series, matched_series_allele_masks, series_match_masks= self.expressed_series_from_allele_mask(allele_mask, min_match)
        return self.match_mask_html(allele_mask, matched_series, matched_series_allele_masks, series_match_masks, 
                                    maybe_allele_mask=allele_mask)

    def match_series_html_from_data_index(self, data_index, min_match=0.0001) :
        allele_mask = self.allele_mask_from_data_index(data_index)
        return self.match_series_html(allele_mask, min_match)
        
    def mask_selected_series_html(self, allele_mask, min_match=0.9, maybe_allele_mask=None) :
        series_data, series_allele_masks, match_allele_masks = self.series_with_allele_mask_matches(allele_mask, min_match)
        return self.match_mask_html(allele_mask, series_data, series_allele_masks, match_allele_masks, maybe_allele_mask)
            