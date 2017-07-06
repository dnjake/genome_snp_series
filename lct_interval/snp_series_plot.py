# -*- coding: utf-8 -*-

import numpy as np
from .series_plot_data import series_plot_data_cls
from .bokeh_series_plot import lct_interval_plot_cls
#from lct_interval.lct_interval_data import lct_data_obj



class superset_yes_no_series_indexes_cls(object) :
    #maybe_allele_mask = lct_data_obj.maybe_allele_mask_from_no_data_indexes(no_data_indexes)
    def __init__(self, plot_context) :
        self.plot_context = plot_context
        
    def do_plot(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        yes_indexes = self.plot_context.yes_data_indexes
        no_indexes = self.plot_context.no_data_indexes
        min_match = self.plot_context.min_allele_match
        #self.focus_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
        exp_data = interval_data_obj.superset_data_from_yes_no_indexes(yes_indexes, no_indexes, min_match)
        self.focus_allele_mask, matched_series, matched_series_allele_masks, series_match_masks = exp_data
        self.plot_context.yes_allele_mask = self.focus_allele_mask
        self.matched_series = matched_series
        self.matched_series_allele_masks = matched_series_allele_masks
        self.series_match_masks = series_match_masks
        for_layout_series_data = interval_data_obj.to_layout_series_data(matched_series)        
        layout_obj = self.plot_context.layout_obj
        layout_obj.set_series_data(for_layout_series_data)
        self.plot_context.layout_data = layout_obj.assign_rows()
        sd = series_plot_data_cls(self.plot_context.chrom, self.plot_context.layout_data)
        self.series_data = sd
        if no_indexes is not None :
            no_allele_mask = interval_data_obj.or_allele_mask(no_indexes)
            self.no_allele_mask = no_allele_mask
            sd.set_no_allele_mask(no_allele_mask)
        sd.generate_plot_data()
        plot_data = sd.plot_data
        self.plot_data = plot_data
        plt = lct_interval_plot_cls(plot_data)
        plt.do_plot()
        return plt.plot_figure

    def get_html(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        maybe_allele_mask = None
        if self.plot_context.no_data_indexes is not None :
            maybe_allele_mask = np.logical_not(self.no_allele_mask)
        html = interval_data_obj.match_mask_html(self.focus_allele_mask, 
                                                 self.matched_series, 
                                                 self.matched_series_allele_masks, 
                                                 self.series_match_masks, 
                                                 maybe_allele_mask=maybe_allele_mask)
        return html                                                 

class subset_yes_no_series_indexes_cls(object) :
    #maybe_allele_mask = lct_data_obj.maybe_allele_mask_from_no_data_indexes(no_data_indexes)
    def __init__(self, plot_context) :
        self.plot_context = plot_context
        
    def do_plot(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        yes_indexes = self.plot_context.yes_data_indexes
        no_indexes = self.plot_context.no_data_indexes
        min_match = self.plot_context.min_allele_match
        self.focus_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
        self.plot_context.yes_allele_mask = self.focus_allele_mask
        #self.focus_allele_mask = match_allele_mask
        exp_data = interval_data_obj.subset_data_from_yes_no_indexes(yes_indexes, no_indexes, min_match)
        matched_series, matched_series_allele_masks, series_match_masks = exp_data
        self.matched_series = matched_series
        self.matched_series_allele_masks = matched_series_allele_masks
        self.series_match_masks = series_match_masks
        for_layout_series_data = interval_data_obj.to_layout_series_data(matched_series)        
        layout_obj = self.plot_context.layout_obj
        layout_obj.set_series_data(for_layout_series_data)
        self.plot_context.layout_data = layout_obj.assign_rows()
        sd = series_plot_data_cls(self.plot_context.chrom, self.plot_context.layout_data)
        self.series_data = sd
        if no_indexes is not None :
            no_allele_mask = interval_data_obj.or_allele_mask(no_indexes)
            self.no_allele_mask = no_allele_mask
            sd.set_no_allele_mask(no_allele_mask)
        sd.generate_plot_data()
        plot_data = sd.plot_data
        self.plot_data = plot_data
        plt = lct_interval_plot_cls(plot_data)
        plt.do_plot()
        return plt.plot_figure
        
    def get_html(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        maybe_allele_mask = None
        if self.plot_context.no_data_indexes is not None :
            maybe_allele_mask = np.logical_not(self.no_allele_mask)
        html = interval_data_obj.match_mask_html(self.focus_allele_mask, 
                                                 self.matched_series, 
                                                 self.matched_series_allele_masks, 
                                                 self.series_match_masks, 
                                                 maybe_allele_mask=maybe_allele_mask)
        return html                                                 
        
class superset_allele_mask_cls(object) :
    #maybe_allele_mask = lct_data_obj.maybe_allele_mask_from_no_data_indexes(no_data_indexes)
    def __init__(self, plot_context) :
        self.plot_context = plot_context
        
    def do_plot(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        yes_allele_mask = self.plot_context.yes_allele_mask
        #yes_indexes = self.plot_context.yes_data_indexes
        #no_indexes = self.plot_context.no_data_indexes
        min_match = self.plot_context.min_allele_match
        #self.focus_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
        self.focus_allele_mask = yes_allele_mask
        #exp_data = interval_data_obj.subset_data_from_yes_no_indexes(yes_indexes,no_indexes, min_match)
        match_data = interval_data_obj.superset_data_from_match_allele_mask(yes_allele_mask, min_match)
        matched_series, matched_series_allele_masks, series_match_masks = match_data
        self.matched_series = matched_series
        self.matched_series_allele_masks = matched_series_allele_masks
        self.series_match_masks = series_match_masks
        for_layout_series_data = interval_data_obj.to_layout_series_data(matched_series)        
        layout_obj = self.plot_context.layout_obj
        layout_obj.set_series_data(for_layout_series_data)
        self.plot_context.layout_data = layout_obj.assign_rows()
        sd = series_plot_data_cls(self.plot_context.chrom, self.plot_context.layout_data)
        self.series_data = sd
        sd.generate_plot_data()
        plot_data = sd.plot_data
        self.plot_data = plot_data
        plt = lct_interval_plot_cls(plot_data)
        plt.do_plot()
        return plt.plot_figure
        
    def get_html(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        html = interval_data_obj.match_mask_html(self.focus_allele_mask, 
                                                 self.matched_series, 
                                                 self.matched_series_allele_masks, 
                                                 self.series_match_masks) 
        return html                                                 




class by_series_indexes_cls(object) :
    def __init__(self, plot_context) :
        self.plot_context = plot_context
    
    def do_plot(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        plot_indexes = self.plot_context.plot_data_indexes
        in_series_data, allele_masks = interval_data_obj.series_data_from_data_indexes(plot_indexes)
        for_layout_series_data = interval_data_obj.to_layout_series_data(in_series_data)        
        layout_obj = self.plot_context.layout_obj
        layout_obj.set_series_data(for_layout_series_data)
        self.plot_context.layout_data = layout_obj.assign_rows()
        sd = series_plot_data_cls(self.plot_context.chrom, self.plot_context.layout_data)
        self.series_data = sd
        if self.plot_context.agg_data_indexes is not None :
            agg_indexes = self.plot_context.agg_data_indexes
            label = self.plot_context.agg_series_label
            allele_count = self.plot_context.agg_series_allele_count
            sd.set_aggregate_snps_series(agg_indexes, label=label, allele_count=allele_count )
        sd.generate_plot_data()
        plot_data = sd.plot_data
        self.plot_data = plot_data
        plt = lct_interval_plot_cls(plot_data)
        plt.do_plot()
        return plt.plot_figure
    
    
class common_series_cls(object) :    
    def __init__(self, plot_context) :
        self.plot_context = plot_context
    
    def do_plot(self) :
        interval_data_obj = self.plot_context.interval_data_obj
        yes_indexes = self.plot_context.yes_data_indexes
        no_indexes = self.plot_context.no_data_indexes
        #min_match = self.plot_context.min_allele_match
        in_series_data, allele_masks = interval_data_obj.series_data_from_data_indexes(yes_indexes)
        for_layout_series_data = interval_data_obj.to_layout_series_data(in_series_data)        
        layout_obj = self.plot_context.layout_obj
        layout_obj.set_series_data(for_layout_series_data)
        self.plot_context.layout_data = layout_obj.assign_rows()
        sd = series_plot_data_cls(self.plot_context.chrom, self.plot_context.layout_data)
        self.series_data = sd
        or_allele_mask = interval_data_obj.or_allele_mask(no_indexes)
        sd.set_no_allele_mask(or_allele_mask)
        seq_allele_mask = interval_data_obj.and_and_not_or_allele_mask(yes_indexes, no_indexes)
        self.common_allele_mask = seq_allele_mask
        sd.set_sequence_series(yes_indexes, seq_allele_mask, 'common_series')
        sd.generate_plot_data()
        plot_data = sd.plot_data
        self.plot_data = plot_data
        plt = lct_interval_plot_cls(plot_data)
        plt.do_plot()
        return plt.plot_figure
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    