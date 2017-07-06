# -*- coding: utf-8 -*-

import numpy as np
from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
from .bokeh_basis_plot import bokeh_basis_plot_cls
from ..series_anal import snp_series as ss
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
cra = country_region_alleles_cls()
anal_context = plot_context_cls()
lct_obj = anal_context.interval_data_obj

'''
I need to be careful about what I want

I need to refresh my memory on the basic operation of the selection

    Each series is associated with a vector of alleles
    
    yes no for data indexes
    
        find vectors for yes series
        
        and them
        
        find vectors for no series
        
        or them
        
        exclued no series from yes series.
        
        That gives the match mask.
        
        The operation is then to match that mask against
        the masks assoicated with all series in the interval
        
            The question is whether the match reaches a threshold
            compared to
            
                The number of chromosomes in the match mask
                (superset)
                
                or
                
                The number of chromosomes in the tested series
                (subset)

What I want to do is incrementally accumulate an allele mask:

    The new processing just has to be at the mask level  I want
    to be able to use an allele mask as input and also to do
    and operations for the no series.
    
    The match operation should not be different.
    But I want to be able to use an allele mask
    as input to the match particularly for superset
    But maybe I also want subset.

One criteria is that the chromosomes don't match all of a set
of series. 
'''
def allele_mask_from_index(data_index) :
    return lct_obj.allele_mask_from_data_index(data_index)
    
def and_allele_mask(data_indexes) :
    return lct_obj.and_allele_mask(data_indexes)
    
def or_allele_mask(data_indexes) :
    return lct_obj.or_allele_mask(data_indexes)

def and_and_not_or_allele_mask(and_data_indexes, not_or_data_indexes) :
    return lct_obj.and_and_not_or_allele_mask(and_data_indexes, not_or_data_indexes)

def not_or_allele_mask(not_or_data_indexes) :
    return lct_obj.not_or_allele_mask(not_or_data_indexes)

def yes_allele_mask_or_indexes(yes_allele_mask, data_indexes):
    return lct_obj.yes_allele_mask_or_indexes(yes_allele_mask, data_indexes)
    
def yes_allele_mask_and_indexes(yes_allele_mask, data_indexes):
    return lct_obj.yes_allele_mask_and_indexes(yes_allele_mask, data_indexes)

def yes_allele_mask_not_or_indexes(yes_allele_mask, data_indexes) :
    return lct_obj.yes_allele_mask_not_or_indexes(yes_allele_mask, data_indexes)

def yes_allele_mask_not_and_indexes(yes_allele_mask, data_indexes) :
    return lct_obj.yes_allele_mask_not_and_indexes(yes_allele_mask, data_indexes)
    
def match_series_html(allele_mask, min_match=0.0001) :
    return lct_obj.match_series_html(allele_mask, min_match)
    #return lct_obj.superset_html_from_mask(allele_mask, min_match)

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
di_9_1023 = 353907

di_31_73 = 353342
di_4_38 = 353450
di_7_20 = 353570
di_6_57 = 353444
di_7_22 = 353315
di_80_38 = 354128
di_7_17 = 353321
di_26_17 = 353258
di_a0_10_17 = 353401
di_a1_10_17 = 353241
di_15_27 = 353459
di_14_43 = 353394
di_6_35 = 353407
di_12_32 = 353252
di_e_4_16 = 353705
di_e_5_16 = 353456
di_5_17 = 353467
di_af_6_19 = 353266
di_s_6_19 = 353568
di_10_25 = 353386
di_5_26 = 353468
di_19_27 = 353280
di_4_27 = 353719
di_7_102 = 353356
di_13_911 = 353807
di_95_176 = 353787
di_51_176 = 354123
di_28_434 = 353614
di_49_136 = 353304
# 74_210
di_5_30 = 353519

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
di_7_90 = 353285
di_29_68 = 353292
di_9_39 = 353559

# eas
di_9_944 = 353372
di_4_1149 = 353849

# afr
di_67_329 = 353331
di_6_32 = 353556
di_180_251 = 353288
di_10_17 = 353241
di_70_166 = 353538
di_219_26 = 353326
di_209_56 = 353276
di_290_16 = 353261
di_22_73 = 353259
di_6_68 = 353247
di_25_27 = 353496
di_7_49 = 353379
di_8_718 = 353312
di_74_210 = 353249
di_22_35 = 353483	
di_9_545 = 353349
di_9_887 = 353729
di_5_588 = 353764
di_147_38 = 353242

def basis_from_allele_mask(pop_code_allele_mask) :
    plt_context = plot_context_cls()
    plt_context.yes_allele_mask = pop_code_allele_mask
    plt_context.yes_data_indexes = None
    plt_context.no_data_indexes = None
    plt_context.basis_obj = plt_context.basis_data_from_allele_mask(plt_context.yes_allele_mask)
    plt_context.basis_series_data = plt_context.basis_obj.basis_series
    plt_context.plot_data_indexes = plt_context.basis_series_data['data_index']
    plt_obj = bokeh_basis_plot_cls(plt_context.interval_data_obj, plt_context.basis_series_data)
    plt_obj.plot_context = plt_context
    return plt_obj    
    
def basis_for_pop(pop_code) :
    pop_code_allele_mask = cra.region_allele_mask(pop_code)
    return basis_from_allele_mask(pop_code_allele_mask)
    
def all_afr() :
    afr_mask = cra.region_allele_mask('afr')
    afx_mask = cra.region_allele_mask('afx')
    all_afr_mask = np.logical_or(afr_mask, afx_mask)
    return basis_from_allele_mask(all_afr_mask)
    
def all_sas() :
    sas_mask = cra.region_allele_mask('sas')
    sax_mask = cra.region_allele_mask('sax')
    all_sas_mask = np.logical_or(sas_mask, sax_mask)
    return basis_from_allele_mask(all_sas_mask)
    
def all_pops() :
    allele_mask = np.empty(5008, '?')
    allele_mask[:] = True
    return basis_from_allele_mask(allele_mask)    
    
    
def superset_from_data_index(data_index, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [data_index]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_yes_no_series_indexes_cls(plt_context)
    return plt_obj

def subset_from_data_index(data_index, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [data_index]
    plt_context.no_data_indexes = None
    plt_context.min_allele_match = min_match
    plt_obj = ssp.subset_yes_no_series_indexes_cls(plt_context)
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


def superset_basis_yes_no(basis_data_index, yes_indexes=None, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    plt_context.yes_data_indexes = [basis_data_index]
    plt_context.get_basis()
    basis_data_item = plt_context.basis_obj.basis_data_from_data_index(basis_data_index)
    yes_allele_mask = basis_data_item['basis_allele_mask']
    if yes_indexes is not None :
        interval_data_obj = plt_context.interval_data_obj
        data_index_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
        yes_allele_mask = np.logical_and(yes_allele_mask, data_index_allele_mask)
    plt_context.basis_data_index = basis_data_index
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.min_allele_match = min_match
    plt_context.yes_allele_mask = yes_allele_mask
    plt_obj = ssp.superset_allele_mask_cls(plt_context)
    return plt_obj

def superset_allele_mask(yes_allele_mask, yes_indexes=None, no_indexes=None, min_match=0.9) :
    plt_context = plot_context_cls()
    if yes_indexes is not None :
        interval_data_obj = plt_context.interval_data_obj
        data_index_allele_mask = interval_data_obj.yes_no_data_indexes(yes_indexes, no_indexes)
        yes_allele_mask = np.logical_and(yes_allele_mask, data_index_allele_mask)
    plt_context.yes_allele_mask = yes_allele_mask
    plt_context.yes_data_indexes = yes_indexes
    plt_context.no_data_indexes = no_indexes
    plt_context.min_allele_match = min_match
    plt_obj = ssp.superset_allele_mask_cls(plt_context)
    return plt_obj

chrom = 2
def series_html(series_index) :
    series_set = ss.snp_series_set(chrom, [series_index])
    series_html = series_set.series_set_html()
    return series_html




