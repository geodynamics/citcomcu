import tables

file = tables.openFile('empty2.h5', mode='w', title='Empty test file')
file.close()
