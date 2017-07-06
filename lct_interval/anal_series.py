# -*- coding: utf-8 -*-

import numpy as np
from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
#from .bokeh_basis_plot import bokeh_basis_plot_cls
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
from ..lct_interval_snp_anal import lct_interval_snp_anal as snp_anal
from . import series_masks as sm

cra = country_region_alleles_cls()
anal_context = plot_context_cls()
lct_obj = anal_context.interval_data_obj

di_13_1696 = 353240
di_9_944 = 353372
di_5_684 = 353498
di_180_251 = 353288
di_123_1561 = 353478
di_10_2206 = 353921
di_64_1575 = 353901
di_7_1868 = 354170
di_7_1760 = 354033
di_32_1361 = 353791
di_5_47 = 353554
di_12_29 = 353243
di_117_1685 = 353244
di_4_815 = 353384
di_4_1699 = 353462
di_6_1503 = 353283
di_26_1414 = 353797
di_81_857 = 353902
di_24_1504 = 354130
di_14_43 = 353394
di_4_815 = 353384
di_193_843 = 353248
di_62_1265 = 353269
di_11_765 = 353380
di_4_911 = 353604
di_39_1014 = 353984
di_10_1218 = 353790
di_13_911 = 353807
di_95_176 = 353787
di_51_176 = 354123
di_74_210 = 353249

# sas
di_8_267 = 353358
di_4_416 = 353511
di_9_1170 = 353906
di_5_1460 = 353814
di_4_1442 = 353925
di_5_976 = 353935
di_4_416 = 353511
di_6_946 = 353504
di_6_713 = 353521

# eas
di_9_944 = 353372
di_4_1149 = 353849

# afr
di_67_329 = 353331
di_23_226 = 353300
di_22_73 = 353259
di_49_136 = 353304
di_28_434 = 353614

# minor 16_1396 selector series
di_18_60 = 353282
di_5_40 = 353536
di_26_17 = 353258
di_5_20 = 353580
di_4_57 = 353987
di_4_34 = 353308


sa_26_1414 = snp_anal.series_anal_cls(di_26_1414, sm.series_data)
sa_117_1685 = snp_anal.series_anal_cls(di_117_1685, sm.series_data)
sa_123_1561 = snp_anal.series_anal_cls(di_123_1561, sm.series_data)
sa_62_1265 = snp_anal.series_anal_cls(di_62_1265, sm.series_data)
sa_193_843 = snp_anal.series_anal_cls(di_193_843, sm.series_data)
sa_13_1696 = snp_anal.series_anal_cls(di_13_1696, sm.series_data)
sa_67_329 = snp_anal.series_anal_cls(di_67_329, sm.series_data)
sa_9_944 = snp_anal.series_anal_cls(di_9_944, sm.series_data)
sa_5_684 = snp_anal.series_anal_cls(di_5_684, sm.series_data)

def superset_yes_no(yes_indexes, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj













































