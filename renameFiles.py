#!/usr/bin/env python

import os

DIRS = ["H:\\Laura\\2018-11-14 700 rad Cs Irradiated Day 4 enteroids\\2 DAY",
        "H:\\Laura\\2018-11-14 700 rad Cs Irradiated Day 4 enteroids\\BY51",
        "H:\\Laura\\2018-11-14 700 rad Cs Irradiated Day 4 enteroids\\BY52",
        "H:\\Laura\\2018-11-14 700 rad Cs Irradiated Day 4 enteroids\\BY53"]

CHANGES = [ ("_c1+2+3+4+5", "_MERGE"),
            ("_c1", "_HOECHST"),
            ("_c2", "_EDU"),
            ("_c3", "_CASP3"),
            ("_c4", "_TUNEL"),
            ("_c5", "_BF") ]

for d in DIRS:
  for oldName in os.listdir(d):
    newName = oldName
    for i in CHANGES:
      newName = newName.replace(*i)
    os.rename(os.path.join(d,oldName), os.path.join(d,newName))
    #print "Renaming %s to %s" % (os.path.join(d,oldName), os.path.join(d,newName))
