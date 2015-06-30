# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 21:51:22 2015

@author: dima
"""

from itertools import groupby, cycle
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_contig_row(dataframe, x = 'Mb', y = 'lod', kind = 'skyscrapers', baselevel = None, 
              colors = ['r','g','b'], msize = 5, fig = None, ax = None, 
              xgrid_step = 5, xgrid_minor_step = 1, len_by_chr = {}, rescale_x = 1,
              normal_font_spacing = 2):
    "uses pieces of code from https://github.com/brentp"
    
    colors = cycle(colors)
    if x not in dataframe.columns:
        x = dataframe.columns[0]
    if y not in dataframe.columns:
        y = dataframe.columns[-1]
    
    round_tick = lambda x_: np.ceil( x_ / xgrid_minor_step) * xgrid_minor_step
    msize = 5
    tmp_x_start_chr = 0
    x_start_chr = {}
    xs = {}
    ys = {}
    cs = {}
    midpoint_by_chr = {}
    _len_by_chr_ = {}
    for seqid, chr_data in dataframe.groupby(level=0):
        if seqid not in len_by_chr:
            _len_by_chr_[seqid] = round_tick( float(chr_data[x].ix[-1]) * rescale_x)
        else:
            _len_by_chr_[seqid] = round_tick( float(len_by_chr[seqid]) * rescale_x)
        x_start_chr[seqid] = round_tick(tmp_x_start_chr)
        region_xs = [float(it) * rescale_x + x_start_chr[seqid] for it in list(chr_data[x]) ]
        xs[seqid] = region_xs
        ys[seqid] = list(chr_data[y])
        cs[seqid] = colors.__next__()
        tmp_x_start_chr += _len_by_chr_[seqid]
        
    chr_order = dataframe.index.unique()
        
    fig = plt.figure()
    ax = fig.add_subplot(111) # ax = fig.add_axes((0.1, 0.09, 0.88, 0.85))
    
    ax.set_ylabel(y)
    
    if baselevel is None:
        ymin = float('Inf')
        for cc in chr_order:
            ymin = min(ymin, min(ys[cc]) )
        baselevel =  np.floor(ymin) - 10
    ymax = -float('Inf')
    for cc in chr_order:
        ymax = max(ymax, max(ys[cc]))
    ymax = np.ceil(ymax)
    
    for cc in chr_order:
        ax.vlines(x_start_chr[cc], baselevel, ymax, colors='k', alpha=0.333)
        if kind == 'skyscrapers' or kind == 'manhattan':
            ax.vlines(xs[cc], baselevel, ys[cc], colors=cs[cc], alpha=0.5)
        if kind == 'connected' or kind == 'line':
            ax.plot(xs[cc], ys[cc], '.-', color=cs[cc], alpha=0.5)
        else:
            ax.scatter(xs[cc], ys[cc], s=msize, c=cs[cc], alpha=0.8, edgecolors='none')
        
    ax.set_xlim(0, xs[chr_order[-1]][-1] )
    ax.set_ylim(baselevel + 10, ymax)
    
    "define the grid"
    
    chr_break_points = np.array(list(x_start_chr.values()))
    chr_break_points.sort()
    
    x_ticks = []
    x_minor_ticks = []
    x_tick_label = []
    first_chr_xtick_index = {}
    chr_xtick_indices = {}
    font_size_coef = {}
    cum_ii = 0
    for cc in chr_order:
        x_minor_grid = np.arange(0, _len_by_chr_[cc], xgrid_minor_step, dtype = int)
        x_grid = np.arange(0, _len_by_chr_[cc], xgrid_step, dtype = int)
        first_chr_xtick_index[cc] = len(x_tick_label)
        x_tick_label.extend([cc] + list(x_grid[1:]))
        x_ticks.extend(cum_ii + x_grid)
        
        x_minor_ticks.extend(cum_ii + x_minor_grid)
        chr_xtick_indices[cc] = first_chr_xtick_index[cc] + np.arange(0, 1+len(x_grid[1:]), dtype = int)
        cum_ii += _len_by_chr_[cc]
        "make the font for small chromosomes tiny"
        if not x_grid[0] and _len_by_chr_[cc] < normal_font_spacing:
            font_size_coef[cc] = _len_by_chr_[cc] / normal_font_spacing
        else:
            font_size_coef[cc] = 1
    "plot the grid"
    ax.xaxis.get_majorticklabels()
    ax.xaxis.grid('on', which='major')
#    ax.xaxis.grid('on', which='minor')
    
    ax.set_xticks(x_ticks, minor = False)
    ax.set_xticklabels(x_tick_label, minor = False)
    
    ax.set_xticks(x_minor_ticks, minor = True)
    
    for cc, xt in first_chr_xtick_index.items():
        plt.setp( ax.xaxis.get_majorticklabels()[xt], rotation=70 )
        if xt>0 and \
          (x_ticks[xt] - x_ticks[xt-1]) < normal_font_spacing and \
          (xt-1) not in first_chr_xtick_index.values():
            x_tick_label[xt-1] = ''
        "make the font for small chromosomes tiny"
        original_size = plt.getp(ax.xaxis.get_majorticklabels()[xt], 'size')
        plt.setp(ax.xaxis.get_majorticklabels()[xt], size = font_size_coef[cc] * original_size )
    for cc in chr_order:
        #ax.xaxis.get_majorticklabels()[ xt ]._color = cs[cc]
        for xx in chr_xtick_indices[cc]:
            ax.xaxis.get_majorticklabels()[ xx ]._color = cs[cc]
    return fig, ax


def plot_contig_col(self,  x = 'Mb', y = ['p_stat', 'p_slct'], x0 = 'x',
             peaks = False, len_by_chr = {}):
#        ax.plot(self['Mb'] , self['p_slct'], 'r.-', label='selection')
#        ax.plot(self['Mb'] , self['p_stat'], 'g-', label='no selection')
    from matplotlib import gridspec


    def find_peaks(y, tolerance = None):
        ddy = np.diff(y, 2)
        signs = np.diff(np.sign(np.diff(y) )) !=0
        if tolerance is None:
            tolerance = np.median(np.abs(ddy[signs])) / 8
        return np.concatenate(([False,], \
        signs * (ddy < -tolerance), [False,]) )
    
    def plot_peaks(x_, y_, xlab_ = None, ax = plt.gca()):
        if type(xlab_) == type(None):
            xlab_ = x_
        "find peaks"
        ind = find_peaks(y_)
        ax.plot(x_[ind] , y_[ind], 'o', markersize = 8, markerfacecolor = 'none',\
        markeredgecolor = 'r')
        for n in range(len(ind)):
            if ind[n]:
                print('peak: %u' % n, file=sys.stderr)
                x, y, xl = x_[n] , y_[n], xlab_[n] 
                ax.text(x, y, "{:,}".format(xl), style='italic')
        print('number of peaks: %u' % sum(ind), file=sys.stderr)
        return
    
    fig = plt.figure() 
    fig.suptitle( (', '.join(y) if type(y) is list else y) + ' for all contigs' )
    
    nsubpl = len(self.contigs)
    nn = [1]
    ax = []
    height = (1-0.2)/ nsubpl
    left = 0.07
    all_widths = 0
    
    for chromosome in self.contigs:
        all_widths += self.chromosome_lengths['Chr1']
    all_widths =  float(max(sr.chromosome_lengths)) / 0.8
    
#        for chromosome in self.contigs:
#            chr_data = self.data.loc[chromosome]        
    for chromosome, chr_data  in self.data.groupby(level=0):
        if not chromosome in self.contigs:
            continue
        width = self.chromosome_lengths[chromosome] / all_widths
        if width < 0.05:
            continue
        
        bottom = (nsubpl - nn[0]) / len(self.contigs)
        ax_box = [left, bottom, width, height]
        ##############################
        ax.append(  fig.add_subplot(nsubpl,1,nn[0]) )
        ax[-1].set_position(ax_box)

        nn[0] += 1
        chr_data.plot(ax = ax[-1], x = x, y = y, 
           label = y)
        if peaks:
            plot_peaks(np.array(chr_data[x]),
                       np.array(chr_data[y[-1] if type([1,2]) is list else y]),
                       np.array(chr_data[x0]), 
                       ax = ax[-1])
        
    ax[-1].set_xlabel('position, Mb', fontsize=12)
    for aa in ax:
        aa.legend('', frameon = False)
    fig.show()
    print('===================================' , file=sys.stderr)
    return fig, ax