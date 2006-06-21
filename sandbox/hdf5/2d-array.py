import Numeric
import tables

def traversal(nx,nz):
    n = 0
    for x in xrange(nx):
        for z in xrange(nz):
            yield n, x, z
            n += 1

file = tables.openFile('2d-array2.h5', mode='w',
                       title='pytables Numeric array!')

box = Numeric.zeros((5,6), typecode=Numeric.Float)

for n,x,z in traversal(5,6):
    box[x,z] = n

file.createArray(file.root, 'box', box, title='Box array!')

file.close()
