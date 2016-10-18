# -*- coding: utf-8 -*-

import numpy as np
total_alleles = 5008

low_alleles = np.zeros(total_alleles, '?')
low_alleles[0::2] = True
high_alleles = np.zeros(total_alleles, '?')
high_alleles[1::2] = True

def sisters_from_low(allele_mask) :
    low_allele_mask = np.logical_and(low_alleles, allele_mask)
    sister_allele_mask = np.zeros(total_alleles, '?')
    sister_allele_mask[1:] = low_allele_mask[:-1]
    return sister_allele_mask

def sisters_from_high(allele_mask) :
    high_allele_mask = np.logical_and(high_alleles, allele_mask)
    sister_allele_mask = np.zeros(total_alleles, '?')
    sister_allele_mask[:-1] = high_allele_mask[1:]
    return sister_allele_mask

def sister_allele_mask(allele_mask) :
    low_sisters = sisters_from_low(allele_mask)
    high_sisters = sisters_from_high(allele_mask)
    sisters = np.logical_or(low_sisters, high_sisters)
    return sisters
    
def have_both_sisters(allele_mask) :
    m_both = allele_mask[::2] == allele_mask[1::2]
    return m_both
    
def allele_mask_for_both_sisters(allele_mask) :
    m_both = have_both_sisters(allele_mask)
    sisters_mask = np.zeros(total_alleles, '?')
    sisters_mask[::2] = m_both
    sisters_mask[1::2] = m_both
    return sisters_mask