# This script needs to run under environment rnatools
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
#cg = forgi.load_rna("examples/input/1y26.fx", allow_many=False)
import forgi.graph.bulge_graph as fgb
import Bio
import forgi.graph.residue as fgr
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv
import math
import matplotlib.pyplot as plt
import numpy as np
import logging
import itertools
import colorsys
from forgi.visual.mplotlib import _annotate_rna_plot, _find_annot_pos_on_circle
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.ndimage import gaussian_filter
def plot_rna_anno(cg, ax=None, offset=(0, 0), text_kwargs={}, backbone_kwargs={},
             basepair_kwargs={}, color=True, lighten=0, annotations={}, color_array = None):
    '''
    Plot an RNA structure given a set of nucleotide coordinates

    .. note::

        This function calls set_axis_off on the axis. You can revert this by
        using ax.set_axis_on() if you like to see the axis.

    :param cg: A forgi.threedee.model.coarse_grain.CoarseGrainRNA structure
    :param ax: A matplotlib plotting area
    :param offset: Offset the plot by these coordinates. If a simple True is passed in, then
                   offset by the current width of the plot
    :param text_kwargs: keyword arguments passed to matplotlib.pyplot.annotate
                        for plotting of the sequence
    :param backbone_kwargs: keyword arguments passed to matplotlib.pyplot.plot
                        for plotting of the backbone links
    :param basepair_kwargs: keyword arguments passed to matplotlib.pyplot.plot
                        for plotting of the basepair links
    :param lighten: Make circles lighter. A percent value where 1 makes
                    everything white and 0 leaves the colors unchanged
    :param annotations: A dictionary {elem_name: string} or None.
                        By default, the element names (e.g. "s0") are plotted
                        next to the element. This dictionary can be used to
                        override the default element names by costum strings.
                        To remove individual annotations, assign an empty string to the key.
                        To remove all annotations, set this to None.

                        .. warning::

                            Annotations are not shown, if there is not enough space.
                            Annotations not shown are logged with level INFO
    :return: (ax, coords) The axes and the coordinates for each nucleotide
    '''

    import RNA
    import matplotlib.colors as mc
    RNA.cvar.rna_plot_type = 1

    coords = []
    #colors = []
    #circles = []

    bp_string = cg.to_dotbracket_string()
    # get the type of element of each nucleotide
    el_string = cg.to_element_string()
    # i.e. eeesssshhhhsssseeee
    el_to_color = {'f': 'orange',
                   't': 'orange',
                   's': 'green',
                   'h': 'blue',
                   'i': 'yellow',
                   'm': 'red'}

    if ax is None:
        ax = plt.gca()

    if offset is None:
        offset = (0, 0)
    elif offset is True:
        offset = (ax.get_xlim()[1], ax.get_ylim()[1])
    else:
        pass

    vrna_coords = RNA.get_xy_coordinates(bp_string)
    # TODO Add option to rotate the plot
    for i, _ in enumerate(bp_string):
        coord = (offset[0] + vrna_coords.get(i).X,
                 offset[1] + vrna_coords.get(i).Y)
        coords.append(coord)
    coords = np.array(coords)
    # First plot backbone
    bkwargs = {"color":"black", "zorder":0}
    bkwargs.update(backbone_kwargs)
    ax.plot(coords[:,0], coords[:,1], **bkwargs)
    # Now plot basepairs
    basepairs = []
    for s in cg.stem_iterator():
        for p1, p2 in cg.stem_bp_iterator(s):
            basepairs.append([coords[p1-1], coords[p2-1]])
    if basepairs:
        basepairs = np.array(basepairs)
        if color:
            c = "red"
        else:
            c = "black"
            bpkwargs = {"color":c, "zorder":0, "linewidth":3}
            bpkwargs.update(basepair_kwargs)
            ax.plot(basepairs[:,:,0].T, basepairs[:,:,1].T, **bpkwargs)
    # Now plot circles
    for i, coord in enumerate(coords):
        if color_array:
            c = color_array[i]
            h,l,s = colorsys.rgb_to_hls(*mc.to_rgb(c))
            if lighten>0:
                l += (1-l)*min(1,lighten)
            else:
                l +=l*max(-1, lighten)
            if l>1 or l<0:
                print(l)
            c=colorsys.hls_to_rgb(h,l,s)
            circle = plt.Circle((coord[0], coord[1]),
                            color=c)
        elif color:
            c = el_to_color[el_string[i]]
            h,l,s = colorsys.rgb_to_hls(*mc.to_rgb(c))
            if lighten>0:
                l += (1-l)*min(1,lighten)
            else:
                l +=l*max(-1, lighten)
            if l>1 or l<0:
                print(l)
            c=colorsys.hls_to_rgb(h,l,s)
            circle = plt.Circle((coord[0], coord[1]),
                            color=c)
        else:
            circle = plt.Circle((coord[0], coord[1]),
                                edgecolor="black", facecolor="white")

        ax.add_artist(circle)
        if cg.seq:
            if "fontweight" not in text_kwargs:
                text_kwargs["fontweight"]="bold"
            ax.annotate(cg.seq[i+1],xy=coord, ha="center", va="center", **text_kwargs )

    all_coords=list(coords)
    ntnum_kwargs = {"color":"gray"}
    ntnum_kwargs.update(text_kwargs)
    for nt in range(10, cg.seq_length, 10):
        # We try different angles
        annot_pos = _find_annot_pos_on_circle(nt, all_coords, cg)
        if annot_pos is not None:
            ax.annotate(str(nt), xy=coords[nt-1], xytext=annot_pos,
                        arrowprops={"width":1, "headwidth":1, "color":"gray"},
                        ha="center", va="center", zorder=0, **ntnum_kwargs)
            all_coords.append(annot_pos)

    _annotate_rna_plot(ax, cg, all_coords, annotations, text_kwargs)
    datalim = ((min(list(coords[:, 0]) + [ax.get_xlim()[0]]),
                min(list(coords[:, 1]) + [ax.get_ylim()[0]])),
               (max(list(coords[:, 0]) + [ax.get_xlim()[1]]),
                max(list(coords[:, 1]) + [ax.get_ylim()[1]])))

    '''
    min_coord = min(datalim[0][0], datalim[0][1])
    max_coord = max(datalim[1][0], datalim[1][1])
    datalim = ((min_coord, min_coord), (max_coord, max_coord))

    print "min_coord:", min_coord
    print "max_coord:", max_coord
    print "datalime:", datalim
    '''

    width = datalim[1][0] - datalim[0][0]
    height = datalim[1][1] - datalim[0][1]

    #ax.set_aspect(width / height)
    ax.set_aspect('equal', 'datalim')
    ax.update_datalim(datalim)
    ax.autoscale_view()
    ax.set_axis_off()

    return (ax, coords)



def read_viennaout(fname):
    ''' read viennaRNA output '''
    with open(fname) as f:
        lines = f.readlines()
    return lines[1].rstrip(), lines[2].split(' ')[0]

if __name__=='__main__':
    seq, struct = read_viennaout(sys.argv[1])
    df = pd.read_csv(sys.argv[2])

    # get structure folded
    bg = fgb.BulgeGraph.from_dotbracket(struct, seq)


    # parameters
    pseudocount = 0.001
    sigma = 1 # gaussian filter
    cmap = cm.YlOrRd


    # calculate relative information
    ipdist = (df['ip1']+df['ip2'])+pseudocount
    ipdist = gaussian_filter(ipdist/ipdist.sum(), sigma = sigma)
    try:
        indist = df['in1']+df['in2']+pseudocount # ENCODE4
    except:
        indist = df['in1']+pseudocount
    indist = gaussian_filter(indist/indist.sum(), sigma = sigma)
    trun = ipdist*np.log(ipdist/indist) ### relative entropy

    f, ax = plt.subplots(figsize = (10,10))
    norm = mpl.colors.Normalize(vmin=trun.min(), vmax=trun.max())
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_anno = [m.to_rgba(x) for x in trun]
    plot_rna_anno(bg, lighten=0.1,
                backbone_kwargs={"linewidth":1}, color = True, 
                color_array = color_anno,
                ax = ax)
    ax.set_title('Relative entropy')

    cb = plt.colorbar(m)
    cb.set_label('Relative entropy')
    plt.savefig(sys.argv[1]+'.png')

