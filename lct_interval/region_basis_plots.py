# -*- coding: utf-8 -*-

import numpy as np
from .lct_interval_plot_context import plot_context_cls
from . import snp_series_plot as ssp
from .bokeh_basis_plot import bokeh_basis_plot_cls
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
cra = country_region_alleles_cls()

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
















    