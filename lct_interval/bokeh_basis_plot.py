# -*- coding: utf-8 -*-

import numpy as np
from bokeh.plotting import figure
from bokeh.models import Range1d, LabelSet, ColumnDataSource
import PyQt4.QtGui as gui
from ..series_anal.snp_series import snp_series_set
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
cra = country_region_alleles_cls()

    
class bokeh_basis_plot_cls(object) :
    #chrom = 2        
    #interval_first_pos = 135757320L
    #interval_last_pos = 136786630L
    gene_dtype = np.dtype([('id', '10S'), ('first_pos', 'i4'), ('last_pos', 'i4'), ('top', 'i4'), ('bottom', 'i4')])
    # starts before interval start('MAP3K19', 135722273, 135782248, 0, 0),
    genes = np.array([ ('RAB3GAP1', 135809835, 135928279, 0, 0), 
                      ('ZRANB3', 135957574, 136288182, 0 ,0),  ('R3HDM1', 136289083, 136482839, 0, 0), 
                      ('UBXN4', 136499189, 136542633, 0, 0), ('LCT', 136545415, 136594750, 0, 0),
                      ('MCM6', 136597196, 136634047, 0, 0), ('DARS', 136664254, 136743222, 0, 0)], gene_dtype)

    series_data_dtype = np.dtype([('data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'), ('row', 'i4'), 
                                 ('label', 'O'), ('allele_count', 'u2'), ('basis_allele_count', 'u2'), ('basis_label', 'O'), ('height', 'i4'), ('top', 'i4'), 
                                 ('bottom', 'i4'), ('color', 'S7'), ('snp_pos', 'O')])
    region_hues = np.array([0, 300, 60, 120, 180 ], dtype='i4')
    gene_height = 70                                 
    bottom_margin = 50
    row_height = 100
    allele_mult = 6
    plot_width = 900
    screen_min_series_len = 20.0
    screen_min_x_pad = 10.0
    screen_max_basis_label_len = 20.0
    y_series_margin = 30
    y_label_margin = 5.0
    toolbar_location = None                                
    def __init__(self, interval_data_obj, basis_data) :
        self.interval_data_obj = interval_data_obj
        self.chrom = interval_data_obj.chrom
        self.interval_first_pos = interval_data_obj.first_pos
        self.interval_last_pos = interval_data_obj.last_pos
        self.interval_series_data = interval_data_obj.series_data
        self.basis_data = basis_data
        self.calc_space_constraints()

    def calc_space_constraints(self) :
        interval_length = float(self.interval_last_pos - self.interval_first_pos)
        sti = interval_length/self.plot_width
        self.min_series_width = int(sti*self.screen_min_series_len)
        self.basis_label_margin = int(sti*self.screen_min_x_pad)
        self.max_basis_label_width = int(sti*self.screen_max_basis_label_len)
        self.max_basis_label_start = self.interval_last_pos - self.max_basis_label_width

    def get_color_for_allele_mask(self, allele_mask) :
        allele_indexes = np.where(allele_mask)[0]
        simple_region_stats = cra.simple_region_stats(allele_indexes)
        obs_to_pred = simple_region_stats['obs_to_pred']
        ind_max_obs_to_pred = obs_to_pred.argmax()
        max_obs_to_pred = obs_to_pred[ind_max_obs_to_pred]
        if max_obs_to_pred > 2.5 :
            sat = 255
            val = 255
        elif max_obs_to_pred > 1.5 :
            sat = 128
            val = 196
        else :
            sat = 0
            val = 128
        qcv = gui.QColor.fromHsv(self.region_hues[ind_max_obs_to_pred], sat, val)
        hrgb = hex(qcv.rgb())
        bokeh_color = '#' + hrgb[4:10]
        return bokeh_color

        series_set_obj = snp_series_set(self.chrom, self.base_series_layouts['data_index'])
        series_set_obj.get_series_and_snp_data()
        series_set_obj.build_struct_data()
        self.base_series_set_data = series_set_obj.set_data

    def build_display_series_data(self) :
        series_data = np.zeros(self.basis_data.size, self.series_data_dtype)
        series_data['data_index'] = self.basis_data['data_index']
        series_data['basis_allele_count'] = self.basis_data['basis_allele_count']
        colors = series_data['color']
        series_allele_masks = self.basis_data['series_allele_mask']
        for ind in range(series_data.size) :
            colors[ind] = self.get_color_for_allele_mask(series_allele_masks[ind])
        self.plot_series_data = series_data
                
    def add_item_data(self, series_data, plot_data) :
        for ind in range(series_data.size) :
            si = series_data[ind]
            pi = plot_data[ind]
            snp_count = si['snp_count']
            allele_count = pi['allele_count']
            pi['label'] = str(snp_count) + '_' + str(allele_count)
            pi['snp_pos'] = si['series_obj'].get_snp_pos()
            pi['basis_label'] = str(pi['basis_allele_count'])

    def add_pos_and_snp_data(self) :
        self.plot_series_data.sort(order='data_index')
        series_set_obj = snp_series_set(self.chrom, self.plot_series_data['data_index'])
        series_set_obj.get_series_and_snp_data()
        series_set_obj.build_struct_data()
        bsd = series_set_obj.set_data
        for name in ['first_pos', 'last_pos', 'allele_count'] :
            self.plot_series_data[name] = bsd[name]
        self.add_item_data(bsd, self.plot_series_data)
       
        
    def set_rows(self) :
        sortargs = np.argsort(self.plot_series_data['allele_count'])
        self.plot_series_data = self.plot_series_data[sortargs]
        row_count = self.plot_series_data.size
        self.plot_series_data['row'] = np.arange(row_count, dtype='i4')
        
    def calc_y_coords(self) :
        heights = np.log2(self.plot_series_data['basis_allele_count']) 
        mask = heights >= 4
        heights[mask] = heights[mask] - 3
        mask = np.logical_not(mask)
        heights[mask] = 1
        heights = self.allele_mult*heights
        heights = heights.astype('i4')
        self.plot_series_data['height'] = heights
        y_coord = self.bottom_margin
        for ind in range(self.plot_series_data.size) :
            ipsd = self.plot_series_data[ind]
            ih = ipsd['height']
            ipsd['bottom'] = y_coord
            ipsd['top'] = y_coord + ih
            y_coord = y_coord + ih + self.y_series_margin
        self.plot_height = y_coord + self.gene_height
        

    '''
    def calc_y_coords(self) :
        self.plot_series_data['bottom'] = self.bottom_margin + self.row_height*self.plot_series_data['row']
        heights = np.log2(self.plot_series_data['basis_allele_count']) 
        mask = heights >= 4
        heights[mask] = heights[mask] - 3
        mask = np.logical_not(mask)
        heights[mask] = 1
        heights = self.allele_mult*heights
        self.plot_series_data['top'] = self.plot_series_data['bottom'] + heights
    '''
    
    def draw_series(self) :
        sd = self.plot_series_data
        m = sd['color'] == '#ffffff'
        if m.any() :
            not_m = np.logical_not(m)
            nsd = sd[not_m]
            sd = sd[m]
            self.plot_figure.quad(left=sd['first_pos'], right=sd['last_pos'] , top=sd['top'], bottom=sd['bottom'], 
                              fill_color=sd['color'], alpha=1.0, line_color='black' )
            sd = nsd
        self.plot_figure.quad(left=sd['first_pos'], right=sd['last_pos'] , top=sd['top'], bottom=sd['bottom'], 
                          fill_color=sd['color'], alpha=1.0, line_color=None)
                
    def draw_series_snps(self) :                              
        lc = 'black'
        colors = []
        x_vals = []
        y_vals = []
        for item in self.plot_series_data :
            item_y_val = (item['bottom'], item['top'])
            item_snp_pos = item['snp_pos']
            for item_pos in item_snp_pos :
                item_pos_x_val = (item_pos, item_pos)
                x_vals.append(item_pos_x_val)
                y_vals.append(item_y_val)
                colors.append(lc)        
        self.plot_figure.multi_line(x_vals, y_vals, color=colors)
                                              
    def draw_series_labels(self) :
        label_x_pos = self.plot_series_data['first_pos']
        label_y_pos = self.plot_series_data['top'] + self.y_label_margin
        label_vals = self.plot_series_data['label']
        label_data = {'x': label_x_pos, 'y': label_y_pos, 'vals': label_vals}
        label_source = ColumnDataSource(label_data)
        labels = LabelSet(x='x', y='y', text='vals', source=label_source, level='glyph', render_mode='canvas',
                          text_baseline='bottom', text_align='left',
                          text_font_size=('9pt'), text_font_style=('bold'), text_alpha=1.0)
        self.plot_figure.add_layout(labels)

    
    def plot_basis_labels(self, x_vals, y_vals, label_vals, text_baseline, alignment) :
        label_data = {'x': x_vals, 'y': y_vals, 'vals': label_vals}
        label_source = ColumnDataSource(label_data)
        labels = LabelSet(x='x', y='y', text='vals', source=label_source, level='glyph', render_mode='canvas',
                          text_baseline=text_baseline, text_align=alignment,
                          text_font_size=('9pt'), text_font_style=('bold'), text_alpha=1.0)
        self.plot_figure.add_layout(labels)

    def draw_basis_labels(self) :
        # make sure there is room at end
        psd = self.plot_series_data
        label_y_pos = psd['bottom'] + ((psd['top'] - psd['bottom'])/2)
        label_x_pos = psd['last_pos'] + self.basis_label_margin
        m_small = label_x_pos - (psd['first_pos']) < self.min_series_width
        label_x_pos[m_small] = psd['first_pos'][m_small] + self.min_series_width
        #max_label_start = self.interval_last_pos - self.min_basis_label_width
        m_need_space = label_x_pos > self.max_basis_label_start
        m_good = np.logical_not(m_need_space)            
        if np.any(m_need_space) :
            label_x_pos[m_need_space] = psd['first_pos'][m_need_space] - self.basis_label_margin
            m_on_top = label_x_pos < self.interval_first_pos
            if np.any(m_on_top) :
                label_x_pos[m_on_top] = psd['last_pos'][m_on_top] - self.max_basis_label_width
                label_y_pos[m_on_top] = psd['top'][m_on_top] + self.y_label_margin
                m_need_space = np.logical_xor(m_need_space, m_on_top)
                self.plot_basis_labels(label_x_pos[m_on_top], label_y_pos[m_on_top],
                                       psd['basis_label'][m_on_top], 'bottom', 'left')
            if np.any(m_need_space) :
                self.plot_basis_labels(label_x_pos[m_need_space], label_y_pos[m_need_space],
                                       psd['basis_label'][m_need_space], 'middle', 'right')
        if np.any(m_good) :            
            self.plot_basis_labels(label_x_pos[m_good], label_y_pos[m_good],
                                   psd['basis_label'][m_good], 'middle', 'left')

    def draw_interval_genes(self) :
        genes = self.genes.copy()
        genes['top'] = self.plot_height - 25
        genes['bottom'] = genes['top'] - 30
        self.plot_figure.quad(left=genes['first_pos'], right=genes['last_pos'], top=genes['top'],
                              bottom=genes['bottom'], fill_color='white', line_color='black')    
        gl_y = 0.5*(genes['top'] - genes['bottom']) + genes['bottom']
        gl_x = 0.5*(genes['last_pos'] - genes['first_pos']) + genes['first_pos']
        gl_source = ColumnDataSource({'x': gl_x, 'y':gl_y, 'vals':genes['id']})
        gene_labels = LabelSet(x='x', y='y', text='vals', source=gl_source, level='glyph', render_mode='canvas',
                               text_baseline='middle', text_align='center', text_font_size=('8pt'), text_alpha=1.0)
        self.plot_figure.add_layout(gene_labels)

    def do_plot(self) :
        self.build_display_series_data()
        self.add_pos_and_snp_data()
        self.set_rows()
        self.calc_y_coords()
        self.plot_figure = figure(plot_width=self.plot_width, plot_height=self.plot_height,
                                  toolbar_location=self.toolbar_location)
        self.plot_figure.yaxis.visible = None
        self.plot_figure.y_range = Range1d(0, self.plot_height)
        self.plot_figure.x_range = Range1d(self.interval_first_pos, self.interval_last_pos)
        self.draw_series()
        self.draw_series_snps()
        self.draw_series_labels()
        self.draw_basis_labels()
        self.draw_interval_genes()
        
