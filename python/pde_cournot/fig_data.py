# SPDX-License-Identifier: MIT
from itertools import cycle

LS = cycle(["solid", "dashed", "dashdot", "dotted",
            (0, (3, 1, 1, 1)),       # densely dashdotted
            (0, (5, 1)),             # densely dashed
            (0, (1, 1)),             # densely dotted
            (0, (3, 1, 1, 1, 1, 1)), # densely dash dotdotted
           ])

color = cycle([None])

# was (16,7.5)
figsize = (12,4.5)
# width/height
figsize_single = (6,4.5)


