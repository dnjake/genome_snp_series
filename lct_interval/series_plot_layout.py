# -*- coding: utf-8 -*-


import numpy as np

        
class by_row_series_layout_cls(object) :
    plot_width_px = 900.0
    min_space_px = 50.0
    work_data_dtype = np.dtype([('data_index', 'u4'), ('first_pos', 'u4'), ('last_pos', 'u4'), 
                                ('allele_count', 'u2'), ('row', 'i4')])
                                    
    out_data_dtype = np.dtype([('data_index', 'u4'), ('row', 'i4')])

    def __init__(self, plot_range) :
        self.plot_first_pos, self.plot_last_pos = plot_range
    def set_series_data(self, series_data) :
        self.work_data = np.zeros(series_data.size, self.work_data_dtype)
        if self.work_data.size == 0 :
            return
        for name in ['data_index', 'first_pos', 'last_pos', 'allele_count'] :
            self.work_data[name] = series_data[name]
        self.work_data.sort(order='allele_count')
        self.work_data = self.work_data[::-1]
        self.find_min_space()
        
    def find_min_space(self) :
        self.min_series_space = self.min_space_px*float(self.plot_last_pos - self.plot_first_pos)/self.plot_width_px                                     
        series_spaces = self.work_data['last_pos'] - self.work_data['first_pos'] + 1
        series_spaces_min = series_spaces.min()
        if series_spaces_min > self.min_series_space:
            self.min_series_space = series_spaces_min
        #print self.min_series_space

    def greater_then_min_spaces(self, in_spaces) :
        out_spaces = []
        for space in in_spaces :
            row, first_pos, last_pos = space
            series_space = last_pos - first_pos + 1
            if series_space > self.min_series_space :
                #print 'series_space', row, first_pos, last_pos, series_space
                out_spaces.append(space)
        if len(out_spaces) > 0 :
            return out_spaces                
        
    def new_spaces(self, old_space, assigned_space) :
        old_row, old_first_pos, old_last_pos = old_space
        new_row, new_first_pos, new_last_pos = assigned_space
        early_space = new_row, old_first_pos, new_first_pos - 1
        late_space = new_row, new_last_pos + 1, old_last_pos
        return self.greater_then_min_spaces((early_space, late_space))

    '''
    would like to assign at least a fixed space for a series except when
    the series is at the end of the plot
    
    also need at least the fixed space for the series except when the series
    is at the end of the plot
    '''

    def assign_series(self, series, spaces) :
        assigned = False
        out_spaces = []
        for space in spaces :
            if assigned :
                out_spaces.append(space)
                continue
            row, space_first_pos, space_last_pos = space
            #print row, space_first_pos, space_last_pos
            if ((series['first_pos'] >= space_first_pos) and 
                (series['last_pos'] <= space_last_pos) and
                ((space_last_pos - series['first_pos'] > self.min_series_space) or
                ((series['first_pos'] + self.min_series_space) > self.plot_last_pos))) :
                assigned = True
                series['row'] = row
                first_pos = series['first_pos']
                last_pos = series['last_pos']
                if last_pos - first_pos < self.min_series_space :
                    last_pos = first_pos + self.min_series_space
                if last_pos > space_last_pos :
                    last_pos = space_last_pos
                    if first_pos > (last_pos - self.min_series_space) :
                        first_pos = last_pos - self.min_series_space
                        if first_pos < space_first_pos :
                            first_pos = space_first_pos
                assigned_space = row, first_pos, last_pos
                #print row, first_pos, last_pos, last_pos - first_pos + 1
                new_spaces = self.new_spaces(space, assigned_space)
                if new_spaces is not None :
                    out_spaces.extend(new_spaces)
        return assigned, out_spaces
        
    def assign_rows(self) :
        row = 0
        row_spaces = [(row, self.plot_first_pos, self.plot_last_pos)]
        need_new_row = False
        for ind in xrange(self.work_data.size) :
            series = self.work_data[ind]
            assigned, out_spaces = self.assign_series(series, row_spaces)
            if assigned :
                if len(out_spaces) > 0 :
                    row_spaces = out_spaces
                else :
                    need_new_row = True
            if (not assigned) or need_new_row :
                row += 1
                row_spaces = [(row, self.plot_first_pos, self.plot_last_pos)]
                need_new_row = False
            if not assigned :
                assigned, out_spaces = self.assign_series(series, row_spaces)
                #print series
                assert assigned
                if len(out_spaces) > 0 :
                    row_spaces = out_spaces
                else :
                    row += 1
                    row_spaces = [(row, self.plot_first_pos, self.plot_last_pos)]
                    need_new_row = False
        out_data = np.zeros(self.work_data.size, self.out_data_dtype)
        out_data['data_index'] = self.work_data['data_index']
        out_data['row'] = self.work_data['row']
        self.out_data = out_data
        return self.out_data
                
                













                                