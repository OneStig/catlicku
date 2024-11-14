# make problem directories a, ..., z
for l in {a..z}; do mkdir $l; cp temp.cpp $l/$l.cpp; done
