# BSD 3-Clause License
# --------------------

# Copyright 2017, The Biotite contributors
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This code is adapted from the biotite package

from sugar.core.fts import Location
import warnings


def _parse_locs(loc_str):
    locs = []
    if loc_str.startswith(("join", "order")):
        warnings.warn('Parsing of genbank loc join, order is untested')
        str_list = loc_str[loc_str.index("(")+1:loc_str.rindex(")")].split(",")
        for s in str_list:
            locs.extend(_parse_locs(s.strip()))
    elif loc_str.startswith("complement"):
        warnings.warn('Parsing of genbank loc complement is untested')
        compl_str = loc_str[loc_str.index("(")+1:loc_str.rindex(")")]
        compl_locs = [
            Location(loc.first, loc.last, Location.Strand.REVERSE, loc.defect)
            for loc in _parse_locs(compl_str)
        ]
        locs.extend(compl_locs)
    else:
        locs = [_parse_single_loc(loc_str)]
    return locs


def _parse_single_loc(loc_str):
    if ".." in loc_str:
        split_char = ".."
        defect = Location.Defect.NONE
    elif "." in loc_str:
        split_char = "."
        defect = Location.Defect.UNK_LOC
    elif "^" in loc_str:
        split_char = "^"
        loc_str_split = loc_str.split("..")
        defect = Location.Defect.BETWEEN
    else:
        # Parse single location
        defect = Location.Defect.NONE
        if loc_str[0] == "<":
            loc_str = loc_str[1:]
            defect |= Location.Defect.BEYOND_LEFT
        elif loc_str[0] == ">":
            loc_str = loc_str[1:]
            defect |= Location.Defect.BEYOND_RIGHT
        first_and_last = int(loc_str)
        return Location(first_and_last-1, first_and_last, defect=defect)
    # Parse location range
    loc_str_split = loc_str.split(split_char)
    first_str = loc_str_split[0]
    last_str = loc_str_split[1]
    # Parse Defects
    if first_str[0] == "<":
        first = int(first_str[1:])
        defect |= Location.Defect.BEYOND_LEFT
    else:
        first = int(first_str)
    if last_str[0] == ">":
        last = int(last_str[1:])
        defect |= Location.Defect.BEYOND_RIGHT
    else:
        last = int(last_str)
    return Location(first-1, last, defect=defect)
