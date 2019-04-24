#!/usr/bin/env python
import sys
import collections
from virtual_family_library import xopen




if __name__=='__main__':
     file_in =sys.argv[1]
     position_dict = collections.defaultdict(list)
     with xopen(file_in,'r') as f:
         for line in f:
             data = line.strip().split("\t")
             if data[0]+data[1] in position_dict:
                 continue
             else:
                 position_dict[data[0]+data[1]].append(line)
                 print line,
    
         
