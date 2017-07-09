# -*- coding: utf-8 -*-
import PyQt4.QtGui as gui
import numpy as np

'''
region_hues = np.array([0, 300, 60, 120, 180 ], dtype='i4')

for v in region_hues :
    qcv = gui.QColor.fromHsv(v, 255, 255)
    print hex(qcv.rgb())
'''
    
pop_colors = '''
<div style="margin-left:100px">
<table>
<tr><td style="background-color:#ff0000;width:50px;"></td><td>African</td></tr>
<tr><td style="background-color:#ff00ff;width:50px;"></td><td>American</td></tr>
<tr><td style="background-color:#ffff00;width:50px;"></td><td>East Asian</td></tr>
<tr><td style="background-color:#00ff00;width:50px;"></td><td>European</td></tr>
<tr><td style="background-color:#00ffff;width:50px;"></td><td>South Asian</td></tr>
</table>
</div>
'''