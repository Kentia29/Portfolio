import sys
from ete3 import Tree

t = Tree(sys.argv[1])
t.set_outgroup(t.get_midpoint_outgroup())
print (t.write())

