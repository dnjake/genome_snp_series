# -*- coding: utf-8 -*-

import numpy as np
from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
from .bokeh_basis_plot import bokeh_basis_plot_cls
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
cra = country_region_alleles_cls()


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

di_31_73 = 353342
di_4_38 = 353450
di_7_20 = 353570
di_6_57 = 353444
di_7_22 = 353315
di_80_38 = 354128
di_7_17 = 353321
di_26_17 = 353258
di_10_17 = 353401
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

di_13_911 = 353807


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











    