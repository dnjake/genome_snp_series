# -*- coding: utf-8 -*-
import numpy as np
from bokeh.plotting import figure
from bokeh.models import Range1d, LabelSet, ColumnDataSource

    
class lct_interval_plot_cls(object) :
    chrom = 2        
    interval_first_pos = 135757320L
    interval_last_pos = 136786630L
    gene_dtype = np.dtype([('id', '10S'), ('first_pos', 'i4'), ('last_pos', 'i4'), ('top', 'i4'), ('bottom', 'i4')])
    # starts before interval start('MAP3K19', 135722273, 135782248, 0, 0),
    genes = np.array([ ('RAB3GAP1', 135809835, 135928279, 0, 0), 
                      ('ZRANB3', 135957574, 136288182, 0 ,0),  ('R3HDM1', 136289083, 136482839, 0, 0), 
                      ('UBXN4', 136499189, 136542633, 0, 0), ('LCT', 136545415, 136594750, 0, 0),
                      ('MCM6', 136597196, 136634047, 0, 0), ('DARS', 136664254, 136743222, 0, 0)], gene_dtype)

    series_data_dtype = np.dtype([('first_pos', 'u4'), ('last_pos', 'u4'), ('row', 'i4'), 
                                 ('label', 'O'), ('allele_count', 'u2'), ('height', 'i4'), ('top', 'i4'), 
                                 ('bottom', 'i4'), ('color', 'S7'), ('snp_pos', 'O')])
    gene_height = 50                                 
    bottom_margin = 50
    row_height = 100
    allele_mult = 6                                
    def __init__(self, plot_data, plot_width=900, toolbar_location=None) :
        self.plot_data_to_series_data(plot_data)
        self.plot_width = plot_width
        self.toolbar_location = toolbar_location
        #self.full_series_data_indexes = None
        #self.comp_source_data_indexes = None
    
    def plot_data_to_series_data(self, plot_data) :
        series_data = np.zeros(plot_data.size, self.series_data_dtype)
        if series_data.size > 0 :
            for name in ['first_pos', 'last_pos', 'row', 'allele_count', 'label', 'color', 'snp_pos'] :
                series_data[name] = plot_data[name]
            max_row = series_data['row'].max()
            series_data['row'] = max_row - series_data['row']
            self.plot_height = (max_row + 1)*self.row_height + self.gene_height + self.bottom_margin
        else :
            self.plot_height = self.gene_height + self.bottom_margin            
        self.series_data = series_data
    
    def calc_y_coords(self) :
        self.series_data['bottom'] = self.bottom_margin + self.row_height*self.series_data['row']
        heights = np.log2(self.series_data['allele_count']) 
        mask = heights >= 4
        heights[mask] = heights[mask] - 3
        mask = np.logical_not(mask)
        heights[mask] = 1
        heights = self.allele_mult*heights
        self.series_data['top'] = self.series_data['bottom'] + heights
    

    def draw_series(self) :
        sd = self.series_data
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
        for item in self.series_data :
            item_y_val = (item['bottom'], item['top'])
            item_snp_pos = item['snp_pos']
            for item_pos in item_snp_pos :
                item_pos_x_val = (item_pos, item_pos)
                x_vals.append(item_pos_x_val)
                y_vals.append(item_y_val)
                colors.append(lc)        
        self.plot_figure.multi_line(x_vals, y_vals, color=colors)
                                              
    def draw_series_labels(self) :
        label_x_pos = self.series_data['first_pos']
        label_y_pos = self.series_data['top'] + 5.0
        label_vals = self.series_data['label']
        label_data = {'x': label_x_pos, 'y': label_y_pos, 'vals': label_vals}
        label_source = ColumnDataSource(label_data)
        labels = LabelSet(x='x', y='y', text='vals', source=label_source, level='glyph', render_mode='canvas',
                          text_baseline='bottom', text_align='left',
                          text_font_size=('9pt'), text_font_style=('bold'), text_alpha=1.0)
        self.plot_figure.add_layout(labels)
                                              
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
        self.calc_y_coords()
        self.plot_figure = figure(plot_width=self.plot_width, plot_height=self.plot_height,
                                  toolbar_location=self.toolbar_location)
        self.plot_figure.yaxis.visible = None
        self.plot_figure.y_range = Range1d(0, self.plot_height)
        self.plot_figure.x_range = Range1d(self.interval_first_pos, self.interval_last_pos)
        self.draw_series()
        self.draw_series_snps()
        self.draw_series_labels()
        self.draw_interval_genes()
        




'''
labels = LabelSet(x=label_x_pos, y=label_y_pos, text=label_vals, render_mode='canvas',
                  text_baseline='bottom', text_align='left', text_font='calibri',
                  text_font_size=('11pt'), text_font_style=('bold'), text_alpha=1.0)
'''                  























    
    
    
    
    
    
    
    
    
    
    
    