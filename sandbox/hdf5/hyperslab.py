import tables

file = tables.openFile('hyperslab.h5', 'r')

# access HDF5 object in file hierarchy
mesh = file.root.proc_mesh

# print representation of proc_mesh
print repr(mesh)
print

print "Accessing z=0 slice from proc_mesh"
print mesh[:,:,0]
print

print "Accessing z=7 slice from proc_mesh"
print mesh[:,:,7]

file.close()
