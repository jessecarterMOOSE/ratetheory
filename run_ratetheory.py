from ratetheory import *

# initialize material
material = Material('Zr')

# add some line sinks
material.add_dislocation_line_sink(1.0e8)
print 'added line sinks, sink dictionary looks like:'
print material.sink_values

# add some small loop sinks
material.add_dislocation_loop_sink(1.0e15, 5.0e-7, name='5_nm_loops')
print 'added 5nm loops, sink dictionary looks like:'
print material.sink_values

# add some bigger loop sinks
material.add_dislocation_loop_sink(2.0e15, 10.0e-7, name='10_nm_loops')
print 'added 10nm loops, sink dictionary looks like:'
print material.sink_values

# retrieve various values
print '5 nm loops sink strength: {:.2e}'.format(material.get_sink_strength('5_nm_loops'))
print '10 nm loops sink strength: {:.2e}'.format(material.get_sink_strength('10_nm_loops'))
print 'sink strength of all loops: {:.2e}'.format(material.get_sink_strength('all_loops'))
print 'all sink strength: {:.2e}'.format(material.get_sink_strength('all'))
print 'all sink strength: {:.2e}'.format(material.get_sink_strength())
print 'all sink strength: {:.2e}'.format(material.sink_strength)
