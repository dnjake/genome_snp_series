# -*- coding: utf-8 -*-
from .lct_interval_plot_context import plot_context_cls
from genomes_dnj.interval_anal.event_stats import chrom_event_stats_cls
from bokeh.plotting import figure
from bokeh.models import Range1d, Label
lct_first_pos = plot_context_cls.interval_data_obj.first_pos
lct_last_pos = plot_context_cls.interval_data_obj.last_pos

def chrom2_stats(field_name) :
    chrom = 2
    chrom_events = chrom_event_stats_cls(chrom)
    stats = chrom_events.stats
    x = stats['event_pos']
    y = stats[field_name]
    y_max = y.max()
    p = figure(plot_width=900, plot_height=600, toolbar_location=None)
    p.y_range = Range1d(0, 1.1*float(y_max))
    p.x_range = Range1d(0, x[-1]+1)
    p.line(x, y)
    lct_inds = stats['event_pos'].searchsorted([lct_first_pos, lct_last_pos])
    lct_stats = stats[lct_inds[0]:lct_inds[1]]
    ind_max = lct_stats[field_name].argmax()
    x_lct_max = float(lct_stats['event_pos'][ind_max])
    y_lct_max = float(lct_stats[field_name][ind_max])
    lct_label = Label(x=x_lct_max, y=y_lct_max, text='LCT', y_offset=10, text_align='center',
                      text_baseline='bottom')
    p.add_layout(lct_label)    
    return p


def chrom2_intvl_stats(field_name, first=135000000, last=138000000) :
    chrom = 2
    chrom_events = chrom_event_stats_cls(chrom)
    stats = chrom_events.stats
    iis = stats['event_pos'].searchsorted([first, last])
    stats = stats[iis[0]:iis[1]]
    x = stats['event_pos']
    y = stats[field_name]
    y_max = y.max()
    p = figure(plot_width=900, plot_height=600, toolbar_location=None)
    p.y_range = Range1d(0, 1.1*float(y_max))
    p.x_range = Range1d(first, last)
    p.line(x, y)
    lct_inds = stats['event_pos'].searchsorted([lct_first_pos, lct_last_pos])
    lct_stats = stats[lct_inds[0]:lct_inds[1]]
    ind_max = lct_stats[field_name].argmax()
    x_lct_max = float(lct_stats['event_pos'][ind_max])
    y_lct_max = float(lct_stats[field_name][ind_max])
    lct_label = Label(x=x_lct_max, y=y_lct_max, text='LCT', y_offset=10, text_align='center',
                      text_baseline='bottom')
    p.add_layout(lct_label)    
    return p





'''
Label(x=70, y=70, x_units='screen' text='Some Stuff', render_mode='css',
      border_line_color='black', border_line_alpha=1.0,
      background_fill_color='white', background_fill_alpha=1.0)
      
citation = Label(x=70, y=70, x_units='screen', y_units='screen',
                 text='Collected by Luke C. 2016-04-01', render_mode='css',
                 border_line_color='black', border_line_alpha=1.0,
                 background_fill_color='white', background_fill_alpha=1.0)
                 
p.add_layout(citation)
'''                 