# -*- coding: utf-8 -*-

from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
di_11_765 = 353380
di_4_911 = 353604
di_6_1503 = 353283
di_4_1699 = 353462
di_26_1414 = 353797
di_64_1575 = 353901
di_10_2206 = 353921
di_7_1868 = 354170
di_7_22 = 353315
di_6_35 = 353407
di_6_57 = 353444
di_5_18 = 353692
di_7_20 = 353570

def lct_agg() :
    plt_context = plot_context_cls()
    plot_data_indexes = [di_11_765, di_4_911, di_6_1503, di_4_1699, di_26_1414, di_64_1575, di_10_2206, di_7_1868]
    plt_context.plot_data_indexes = plot_data_indexes
    plt_context.agg_data_indexes = plot_data_indexes
    plt_context.agg_series_label = 'lct_snps'
    plt_context.agg_series_allele_count = 765
    plt_obj = ssp.by_series_indexes_cls(plt_context)
    return plt_obj

def subset_4_1699(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_4_1699] 
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_6_1503(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_6_1503]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_6_1503_26_1414(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_6_1503, di_26_1414]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_4_911(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_4_911]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_6_1503_26_1414_not_4_911(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_6_1503, di_26_1414]
    plt_context.no_data_indexes = [di_4_911, di_6_57, di_6_35, di_7_22, di_7_20, di_11_765]
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def superset_6_1503_26_1414_not_4_911(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_6_1503, di_26_1414]
    plt_context.no_data_indexes = [di_4_911, di_6_57, di_6_35, di_7_22, di_7_20, di_11_765]
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def superset_4_911_not_11_765(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_4_911]
    plt_context.no_data_indexes = [di_11_765]
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def superset_26_1414(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_26_1414]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def superset_11_765(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_11_765]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_11_765(min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [di_11_765]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
    return plt_obj
