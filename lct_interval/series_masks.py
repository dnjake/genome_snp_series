# -*- coding: utf-8 -*-

import numpy as np
from . import series_plots as dm 
from ..lct_interval_snp_anal.lct_interval_snp_anal import snp_and_series_cls
series_data = snp_and_series_cls()

am_117_1685 = series_data.series_allele_mask_from_index(dm.di_117_1685)
am_123_1561 = series_data.series_allele_mask_from_index(dm.di_123_1561)
am_62_1265 = series_data.series_allele_mask_from_index(dm.di_62_1265)
nam_117_1685 = np.logical_not(am_117_1685)
nam_123_1561 = np.logical_not(am_123_1561)
nam_62_1265 = np.logical_not(am_62_1265)
am_193_843 = series_data.series_allele_mask_from_index(dm.di_193_843)
nam_193_843 = np.logical_not(am_193_843)
am_4_815 = series_data.series_allele_mask_from_index(dm.di_4_815)
am_6_946 = series_data.series_allele_mask_from_index(dm.di_6_946)
am_13_1696 = series_data.series_allele_mask_from_index(dm.di_13_1696)
nam_13_1696 = np.logical_not(am_13_1696)
am_117_1685_or_123_1561 = np.logical_or(am_117_1685, am_123_1561)
am_117_1685_or_62_1265 = np.logical_or(am_117_1685, am_62_1265)
am_117_1685_or_123_1561_or_62_1265 = np.logical_or(am_117_1685_or_123_1561, am_117_1685_or_62_1265)
nam_117_1685_or_123_1561_or_62_1265 = np.logical_not(am_117_1685_or_123_1561_or_62_1265)
am_117_1685_not_123_1561 = np.logical_and(am_117_1685, nam_123_1561)
am_123_1561_not_117_1685 = np.logical_and(am_123_1561, nam_117_1685)
am_117_1685_not_62_1265 = np.logical_and(am_117_1685, nam_62_1265)
am_117_1685_not_193_843 = np.logical_and(am_117_1685, nam_193_843)
am_117_1685_not_123_1561_or_13_1696 = np.logical_and(am_117_1685_not_123_1561, nam_13_1696)
am_117_1685_or_123_1561_not_62_1265 = np.logical_and(am_117_1685_or_123_1561, nam_62_1265)
am_4_815_or_all = np.logical_or(am_4_815, am_117_1685_or_123_1561_or_62_1265)
nam_4_815_or_all = np.logical_not(am_4_815_or_all)
am_9_944 = series_data.series_allele_mask_from_index(dm.di_9_944)
am_5_684 = series_data.series_allele_mask_from_index(dm.di_5_684)
am_67_329 = series_data.series_allele_mask_from_index(dm.di_67_329)
am_180_251 = series_data.series_allele_mask_from_index(dm.di_180_251)
am_70_166 = series_data.series_allele_mask_from_index(dm.di_70_166)
am_219_26 = series_data.series_allele_mask_from_index(dm.di_219_26)
am_290_16 = series_data.series_allele_mask_from_index(dm.di_290_16)
am_22_73 = series_data.series_allele_mask_from_index(dm.di_22_73)
am_209_56 = series_data.series_allele_mask_from_index(dm.di_209_56)
am_74_210 = series_data.series_allele_mask_from_index(dm.di_74_210)
am_28_434 = series_data.series_allele_mask_from_index(dm.di_28_434)
am_28_434_117_1685 = np.logical_and(am_28_434, am_117_1685)
am_28_434_117_1685_123_1561 = np.logical_and(am_28_434_117_1685, am_123_1561)
am_28_434_117_1685_not_123_1561 = np.logical_and(am_28_434, am_117_1685_not_123_1561)
am_28_434_117_1685_not_123_1561_13_1696 = np.logical_and(am_28_434_117_1685_not_123_1561, nam_13_1696)
am_28_434_117_1685_13_1696 = np.logical_and(am_28_434_117_1685, am_13_1696)
am_28_434_123_1561 = np.logical_and(am_28_434, am_123_1561)
am_28_434_123_1561_not_117_1685 = np.logical_and(am_28_434_123_1561, nam_117_1685)
am_28_434_123_1561_13_1696_not_117_1685 = np.logical_and(am_28_434_123_1561_not_117_1685, am_13_1696)
nam_9_944 = np.logical_not(am_9_944)
nam_5_684 = np.logical_not(am_5_684)
nam_67_329 = np.logical_not(am_67_329)
nam_180_251 = np.logical_not(am_180_251)
nam_70_166 = np.logical_not(am_70_166)
nam_219_26 = np.logical_not(am_219_26)
nam_290_16 = np.logical_not(am_290_16)
nam_22_73 = np.logical_not(am_22_73)
nam_209_56 = np.logical_not(am_209_56)










