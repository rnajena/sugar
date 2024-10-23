# (C) 2024, Tom Eulenfeld, MIT license

from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Ellipse
from matplotlib.transforms import Bbox
import random
from sugar import BioSeq
from sugar.data import gcode


def create_logo(fname, seed=None, **kw):
    name = 'SUGAR'

    random.seed(seed or datetime.today().toordinal())
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
                      edgecolor='0.7', fc='None', lw=2)
    arc = Arc(xy=(0.5, 0), width=3, height=0.75, theta1=40, theta2=140,
              color='0.7', lw=2)
    ax.add_patch(ellipse)
    ax.add_patch(arc)
    ax.annotate(aa.replace('*', '    '), (0.5, 0.5), ha='center', size=50)
    ax.annotate(seq[i:j], (0.5, 0.44), ha='center', size=32)
    ax.annotate(seq[:i], (0.37, 0.4), ha='right', size=18, color='0.7')
    ax.annotate(seq[j:], (0.63, 0.4), ha='left', size=18, color='0.7')
    fig.savefig(fname, bbox_inches=Bbox([[0, 0.5], [4, 2.5]]), **kw)


if __name__ == '__main__':
    create_logo('sugar_logo.png')
    create_logo('sugar_logo_transparent.png', transparent=True)
