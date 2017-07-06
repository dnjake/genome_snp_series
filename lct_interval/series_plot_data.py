# -*- coding: utf-8 -*-

import numpy as np
from ..series_anal.snp_series import snp_series_set
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
#import html_display.array_table as html
cra = country_region_alleles_cls()

import PyQt4.QtGui as gui

class series_plot_data_cls(object) :
    plot_data_dtype = np.dtype([('first_pos', 'u4'), ('last_pos', 'u4'),  ('snp_count', 'u2'),
                                   ('allele_count', 'u2'), ('row', 'u2'), ('label', 'O'), ('color', 'S7'), 
                                   ('snp_pos', 'O')])
    allele_masks_dtype = np.dtype([('data_index', 'u4'), ('allele_mask', '?', 5008)])
    series_layout_dtype = np.dtype([('data_index', 'u4'), ('row', 'i4')])                                   
    region_hues = np.array([0, 300, 60, 120, 180 ], dtype='i4')
    #chrom = 2            
    def __init__(self, chrom, series_layout) :
        self.chrom = chrom
        self.base_series_layouts = series_layout
        self.base_series_layouts.sort(order='data_index')
        series_set_obj = snp_series_set(self.chrom, self.base_series_layouts['data_index'])
        series_set_obj.get_series_and_snp_data()
        series_set_obj.build_struct_data()
        self.base_series_set_data = series_set_obj.set_data
        self.has_yes_allele_mask = False
        self.has_no_allele_mask = False
        self.has_sequence_series = False
        self.has_aggregate_snps_series = False
        
    def build_base_series_allele_masks(self) :
        allele_masks = []
        for item in self.base_series_set_data :
            obj = item['series_obj']
            item_data = (obj.data_index, obj.series_allele_mask)
            allele_masks.append(item_data)
        self.base_series_allele_masks = np.array(allele_masks, self.allele_masks_dtype)
        
    def set_yes_allele_mask(self, allele_mask) :
        self.yes_allele_mask = allele_mask
        self.has_yes_allele_mask = True
        
    def set_no_allele_mask(self, allele_mask) :
        self.no_allele_mask = allele_mask
        self.has_no_allele_mask = True
        self.build_regions_after_no_allele_mask()
        
    def build_regions_after_no_allele_mask(self) :
        regions = cra.regions.copy()
        maybe_allele_mask = np.logical_not(self.no_allele_mask)
        region_allele_masks = regions['allele_mask']
        for ind in range(region_allele_masks.shape[0]) :
            region_allele_masks[ind] = np.logical_and(region_allele_masks[ind], maybe_allele_mask)
        maybe_alleles_per_region = regions['allele_mask'].sum(axis=1)
        maybe_alleles_total = maybe_alleles_per_region.sum()
        self.maybe_regions = regions
        self.maybe_alleles_per_region = maybe_alleles_per_region.astype('f4')
        self.maybe_alleles_total = float(maybe_alleles_total)
        self.freq_maybe_alleles_per_region = self.maybe_alleles_per_region/self.maybe_alleles_total
        self.maybe_simple_alleles_per_region = np.zeros(cra.simple_region_codes.size, dtype='f4')
        self.maybe_simple_alleles_per_region[[1, 2, 3]] = self.maybe_alleles_per_region[[2, 3, 4]]
        self.maybe_simple_alleles_per_region[0] = self.maybe_alleles_per_region[0] + self.maybe_alleles_per_region[1]
        self.maybe_simple_alleles_per_region[4] = self.maybe_alleles_per_region[5] + self.maybe_alleles_per_region[6]
        self.freq_maybe_simple_alleles_per_region = self.maybe_simple_alleles_per_region/self.maybe_alleles_total
        

    def maybe_region_stats(self, allele_indexes) :
        sample_region_alleles = self.maybe_regions['allele_mask'][:, allele_indexes]
        results = np.zeros(cra.regions.size, cra.analysis_result_dtype)
        results['code'] = cra.region_codes
        results['allele_count'] = sample_region_alleles.sum(axis=1)
        results['freq'] = results['allele_count']/self.maybe_alleles_per_region
        pred = float(allele_indexes.size)*self.freq_maybe_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results

    def maybe_simple_region_stats(self, allele_indexes) :
        trs = self.maybe_region_stats(allele_indexes)
        results = np.zeros(cra.simple_region_codes.size, cra.analysis_result_dtype)
        results['code'] = cra.simple_region_codes
        trs_counts = trs['allele_count']
        results_counts = results['allele_count']
        results_counts[[1, 2, 3]] = trs_counts[[2, 3, 4]]
        results_counts[0] = trs_counts[0] + trs_counts[1]
        results_counts[4] = trs_counts[5] + trs_counts[6]
        results['freq'] = results['allele_count']/self.maybe_simple_alleles_per_region
        pred = float(allele_indexes.size)*self.freq_maybe_simple_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results
        

    def build_display_allele_masks(self) :
        self.display_allele_masks = self.base_series_allele_masks.copy()
        allele_masks = self.display_allele_masks['allele_mask']
        if self.has_yes_allele_mask :
            allele_masks[:] = np.logical_and(allele_masks, self.yes_allele_mask)
        '''
        if self.has_no_allele_mask :
            good_allele_mask = np.logical_not(self.no_allele_mask)            
            allele_masks[:] = np.logical_and(allele_masks, good_allele_mask)
            self.build_regions_after_no_allele_mask()
        '''
         
    def to_nd_array(self, data_indexes) :
        if type(data_indexes) is np.ndarray :
            return data_indexes
        else :
            return np.array(data_indexes, dtype='u4')

    '''
    In the sequence case, I want the colors to be whatever they would be without the sequence
    I want to use the sequence data indexes to figure out the first and last positions for sequence
    series and the snps to draw in it.  The color for the series should come from its allele mask
    in the standard way
    '''

    def set_sequence_series(self, sequence_data_indexes, sequence_allele_mask, 
                            series_display_label='') :
        self.sequence_series_data_indexes = self.to_nd_array(sequence_data_indexes)
        self.sequence_series_allele_mask = sequence_allele_mask
        #self.sequence_series_display_row = series_display_row
        self.sequence_series_display_label = series_display_label
        self.has_sequence_series = True
        
    def set_aggregate_snps_series(self, source_data_indexes, display_color='#ffffff',
                                  label='', allele_count=None) :
        self.aggregate_snps_series_data_indexes = self.to_nd_array(source_data_indexes)
        #self.aggregate_snps_series_display_row = series_display_row
        self.aggregate_snps_series_display_color = display_color
        self.aggregate_snps_series_label = label
        self.aggregate_snps_series_allele_count = allele_count
        self.has_aggregate_snps_series = True
            
    def count_plot_items(self) :
        count = self.base_series_layouts.size
        if self.has_sequence_series :
            count += 1
        if self.has_aggregate_snps_series :
            count += 1
        return count

    def get_source_series_data(self, plot_data, source_data_indexes) :
        source_data_indexes.sort()
        inds = self.base_series_layouts['data_index'].searchsorted(source_data_indexes)
        source_data = plot_data[inds]
        #self.source_data = source_data
        first_pos = source_data['first_pos'].min()
        last_pos = source_data['last_pos'].max()
        snp_pos = []
        for snps in source_data['snp_pos'] :
            snp_pos.append(snps)
        snp_pos = np.concatenate(snp_pos)
        return first_pos, last_pos, snp_pos

    def get_color_for_allele_mask(self, allele_mask) :
        allele_indexes = np.where(allele_mask)[0]
        '''
        I want to stick with a single color
        if self.has_no_allele_mask :
            simple_region_stats = self.maybe_simple_region_stats(allele_indexes)
            #print simple_region_stats
            #print self.maybe_alleles_per_region
            #print self.maybe_alleles_total
        else :
            simple_region_stats = cra.simple_region_stats(allele_indexes)
        '''
        simple_region_stats = cra.simple_region_stats(allele_indexes)        
        obs_to_pred = simple_region_stats['obs_to_pred']
        ind_max_obs_to_pred = obs_to_pred.argmax()
        max_obs_to_pred = obs_to_pred[ind_max_obs_to_pred]
        if max_obs_to_pred > 2.5 :
            sat = 255
            val = 255
        elif max_obs_to_pred > 1.5 :
            sat = 128
            val = 196
        else :
            sat = 0
            val = 128
        qcv = gui.QColor.fromHsv(self.region_hues[ind_max_obs_to_pred], sat, val)
        hrgb = hex(qcv.rgb())
        bokeh_color = '#' + hrgb[4:10]
        return bokeh_color
        
    def plot_data_item(self, layout, series_set_data, display_allele_mask, base_allele_mask) :
        first_pos = series_set_data['first_pos']
        last_pos = series_set_data['last_pos']
        snp_count = series_set_data['snp_count']
        allele_count = display_allele_mask.sum()
        base_allele_count = base_allele_mask.sum()
        row = layout['row']
        label = str(snp_count) + '_' + str(base_allele_count)
        '''
        I want to stick with a single color for a series
        color = self.get_color_for_allele_mask(display_allele_mask)
        '''
        color = self.get_color_for_allele_mask(base_allele_mask)
        snp_pos = series_set_data['series_obj'].get_snp_pos()
        return first_pos, last_pos, snp_count, allele_count, row, label, color, snp_pos

    def get_base_series_plot_data(self, plot_data) :
        for ind in range(self.base_series_layouts.size) :
            layout_item = self.base_series_layouts[ind]
            series_set_data_item = self.base_series_set_data[ind]
            allele_mask_item = self.display_allele_masks[ind]
            base_allele_mask_item = self.base_series_allele_masks[ind]
            plot_data[ind] = self.plot_data_item(layout_item, series_set_data_item, 
                            allele_mask_item['allele_mask'], base_allele_mask_item['allele_mask'])

    def get_aggregate_snp_series_data(self, base_plot_data) :
        source_data_indexes = self.aggregate_snps_series_data_indexes
        first_pos, last_pos, snp_pos = self.get_source_series_data(base_plot_data, source_data_indexes)
        #self.snp_pos = snp_pos
        base = self.base_series_layouts
        base_max_row = base['row'].max()
        row = base_max_row + 1
        #row = self.aggregate_snps_series_display_row
        color = self.aggregate_snps_series_display_color
        label = self.aggregate_snps_series_label
        allele_count = self.aggregate_snps_series_allele_count
        snp_count = snp_pos.size
        return first_pos, last_pos, snp_count, allele_count, row, label, color, snp_pos
        
    def get_sequence_series_data(self, base_plot_data) :
        source_data_indexes = self.sequence_series_data_indexes
        first_pos, last_pos, snp_pos = self.get_source_series_data(base_plot_data, source_data_indexes)
        #self.snp_pos = snp_pos
        base = self.base_series_layouts
        base_max_row = base['row'].max()
        row = base_max_row + 1
        #row = self.sequence_series_display_row
        #print self.sequence_series_allele_mask.sum()        
        color = self.get_color_for_allele_mask(self.sequence_series_allele_mask)
        label = self.sequence_series_display_label
        allele_count = self.sequence_series_allele_mask.sum()
        snp_count = snp_pos.size
        return first_pos, last_pos, snp_count, allele_count, row, label, color, snp_pos

    def generate_plot_data(self) :
        self.build_base_series_allele_masks()
        self.build_display_allele_masks()
        plot_data_size = self.count_plot_items()
        plot_data = np.zeros(plot_data_size, self.plot_data_dtype)
        self.plot_data = plot_data
        self.get_base_series_plot_data(plot_data)
        base_data_size = self.base_series_layouts.size
        next_ind = -1
        base_plot_data = plot_data[:base_data_size]
        if self.has_aggregate_snps_series :
            plot_data[next_ind] = self.get_aggregate_snp_series_data(base_plot_data)
            next_ind -= 1            
        if self.has_sequence_series :
            plot_data[next_ind] = self.get_sequence_series_data(base_plot_data)


