# (C) 2024, Tom Eulenfeld, MIT license
"""
Create the sugar logo
"""
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Ellipse
from matplotlib.transforms import Bbox
import random
from sugar import BioSeq
from sugar.data import gcode


def create_logo(fname, seed=None, color1='k', color2='0.7', **kw):
    name = 'SUGAR'
    try:
        seed = datetime.fromisoformat(seed).toordinal()
    except Exception:
        if seed is None:
            seed = datetime.today().toordinal()
    random.seed(seed)
    gc = gcode()
    first_codon = random.choice(sorted(gc.ttinv[name[0]]))
    last_codon = random.choice(sorted(gc.ttinv[name[-1]]))
    seq = BioSeq(first_codon + name[1:-1] + last_codon).str.replace('T', 'U')
    i, j = seq.match('stop').span()
    aa = seq.copy().translate(complete=True)

    fig = plt.figure(figsize=(4,3), frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.spines[:].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False)
    ellipse = Ellipse(xy=(0.5, 0.5), width=0.8, height=0.5,
                      edgecolor=color2, fc='None', lw=2)
    arc = Arc(xy=(0.5, 0), width=3, height=0.75, theta1=40, theta2=140,
              color=color2, lw=2)
    ax.add_patch(ellipse)
    ax.add_patch(arc)
    ax.annotate(aa.str.replace('*', '    '), (0.5, 0.5), ha='center', size=50, color=color1)
    ax.annotate(seq[i:j], (0.5, 0.44), ha='center', size=32, color=color1)
    ax.annotate(seq[:i], (0.37, 0.4), ha='right', size=18, color=color2)
    ax.annotate(seq[j:], (0.63, 0.4), ha='left', size=18, color=color2)
    if fname is None:
        plt.show()
    else:
        fig.savefig(fname, bbox_inches=Bbox([[0, 0.5], [4, 2.5]]), **kw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fname', nargs='?', help='output filename')
    parser.add_argument('--seed', help='random seed (may be date in isoformat, default is today)')
    parser.add_argument('--color1', help='first color', default='k')
    parser.add_argument('--color2', help='second color', default='0.7')
    parser.add_argument('--transparent', action='store_true', help='use transparent background')
    args = vars(parser.parse_args())
    create_logo(**args)


if __name__ == '__main__':
    main()
