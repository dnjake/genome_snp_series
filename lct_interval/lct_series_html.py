# -*- coding: utf-8 -*-

import genomes_dnj.series_anal.snp_series as ss

chrom = 2

di_11_765 = 353380
di_4_911 = 353604
di_6_1503 = 353283
di_4_1699 = 353462
di_26_1414 = 353797
di_10_2206 = 353921
di_64_1575 = 353901
di_7_1868 = 354170

data_indexes = [di_11_765, di_4_911, di_6_1503, di_4_1699, di_26_1414, di_10_2206, di_64_1575, di_7_1868]

series_set = ss.snp_series_set(chrom, data_indexes)
lct_html = series_set.series_set_html()
