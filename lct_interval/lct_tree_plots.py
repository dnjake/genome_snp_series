# -*- coding: utf-8 -*-

from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
from .bokeh_basis_plot import bokeh_basis_plot_cls

di_26_1414 = 353797
di_117_1685 = 353244
di_123_1561 = 353478
di_62_1295 = 353269
di_209_56 = 353276
di_4_1699 = 353462
di_6_1503 = 353283
di_64_1575 = 353901
di_12_32 = 353252
di_5_16 = 353821
di_10_43 = 353347
di_14_48 = 353833
di_80_38 = 354128
di_4_1149 = 353849
di_4_1442 = 353925
di_7_1760 = 354033
di_24_1504 = 354130
di_32_1361 = 353791
di_6_35 = 353407
di_4_911 = 353604
di_11_765 = 353380
di_10_2206 = 353921
di_7_1868 = 354170
di_13_1696 = 353240
di_95_176 = 353787
di_51_176 = 354123
di_5_18 = 353692

def subset_4_1699(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_4_1699]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def basis_4_1699() :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_4_1699]
    plt_context.no_data_indexes = None
    plt_context.basis_obj = plt_context.context_basis_data()
    plt_context.basis_series_data = plt_context.basis_obj.basis_series
    plt_context.plot_data_indexes = plt_context.basis_series_data['data_index']
    plt_obj = bokeh_basis_plot_cls(plt_context.interval_data_obj, plt_context.basis_series_data)
    plt_obj.plot_context = plt_context
    return plt_obj

def superset_yes_no(yes_indexes, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj


def subset_yes_no(yes_indexes, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def basis_yes_no(yes_indexes, no_indexes=None) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.basis_obj = plt_context.context_basis_data()
    plt_context.basis_series_data = plt_context.basis_obj.basis_series
    plt_context.plot_data_indexes = plt_context.basis_series_data['data_index']
    #plt_obj = bokeh_basis_plot_cls(plt_context.interval_data_obj, plt_context.basis_series_data)
    #plt_obj.plot_context = plt_context
    return plt_context




















