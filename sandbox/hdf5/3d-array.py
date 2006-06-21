import Numeric
import tables

def traversal(nx,ny,nz):
    n = 0
    for y in xrange(ny):
        for x in xrange(nx):
            for z in xrange(nz):
                yield n, x, y, z
                n += 1

file = tables.openFile('3d-array2.h5', mode='w', 
                       title='pytables Numeric array!')

box = Numeric.zeros((5,6,2), typecode=Numeric.Float)

for n,x,y,z in traversal(5,6,2):
    box[x,y,z] = n

file.createArray(file.root, 'box', box, title='Box array!')
file.close()
