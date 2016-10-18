# -*- coding: utf-8 -*-

from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
from ..series_anal import by_haplotype as sis
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
di_11_765 = 353380
di_4_911 = 353604
di_10_2206 = 353921
di_7_1868 = 354170


def sister_yes_no(yes_indexes, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    interval_data_obj = plt_context.interval_data_obj
    plt_context.initial_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
    sister_allele_mask = sis.sister_allele_mask(plt_context.initial_allele_mask)
    #plt_context.yes_data_indexes = yes_indexes
    #plt_context.no_data_indexes = no_indexes
    plt_context.yes_allele_mask = sister_allele_mask
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_allele_mask_cls(plt_context)
    return plt_obj
