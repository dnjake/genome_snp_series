# -*- coding: utf-8 -*-

from ..series_anal.series_relation_data import interval_data_cls
from .series_plot_layout import by_row_series_layout_cls
from ..series_anal.basis_series import basis_series_cls
from .bokeh_basis_plot import bokeh_basis_plot_cls

chrom = 2        
interval_first_pos = 135757320L
interval_last_pos = 136786630L
lct_data_obj = interval_data_cls(chrom, interval_first_pos, interval_last_pos)

class plot_context_cls(object) :
    chrom = 2        
    interval_first_pos = 135757320L
    interval_last_pos = 136786630L
    interval_mid_pos = (interval_last_pos - interval_first_pos)/2
    interval_data_obj = interval_data_cls(chrom, interval_first_pos, interval_last_pos)
    plot_range = interval_data_obj.first_pos, interval_data_obj.last_pos
    yes_data_indexes = None
    no_data_indexes = None
    agg_data_indexes = None
    basis_obj = None
    yes_allele_mask = None

    def __init__(self, layout_obj=None) :
        self.layout_obj = layout_obj
        if self.layout_obj is None :
            self.layout_obj = by_row_series_layout_cls(self.plot_range)

    def basis_data_from_allele_mask(self, allele_mask) :
        bs_obj = basis_series_cls(self.interval_data_obj, allele_mask)
        return bs_obj
        
    def context_basis_data(self) : 
        basis_allele_mask = self.interval_data_obj.yes_no_data_indexes(self.yes_data_indexes, self.no_data_indexes)
        return self.basis_data_from_allele_mask(basis_allele_mask)
        
    def get_basis(self) :
        self.basis_obj = self.context_basis_data()
        self.basis_series_data = self.basis_obj.basis_series
        self.basis_plot_data_indexes = self.basis_series_data['data_index']
        
    def get_basis_plot(self) :
        self.get_basis()        
        plt_obj = bokeh_basis_plot_cls(self.interval_data_obj, self.basis_series_data)
        plt_obj.plot_context = self
        return plt_obj

    def do_basis_plot(self) :
        basis_plt_obj = self.get_basis_plot()
        basis_plt_obj.do_plot()
        basis_plt_obj.plot_context = None
        return basis_plt_obj.plot_figure
        
    def get_basis_html(self) :
        if self.basis_obj is None :
            self.get_basis()
        return self.basis_obj.html_sorted_by_allele_count()
        
    def get_country_html(self) :
        if self.yes_allele_mask is not None :
            return self.interval_data_obj.cra.country_html_table(self.yes_allele_mask)
        
        
        
        
        
        
        
        
        
        
        
        
        
        