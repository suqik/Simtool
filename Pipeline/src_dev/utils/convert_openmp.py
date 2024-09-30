""" This script converts a FastPM snapshot to a Gadget-1 snapshot.

    Three blocks are converted, position, velocity and ID>

    precision of position and velocity is controled via --precision
    ID is always 8 bytes long.

    The default is 1 million particles per file. Change it to something bigger
"""

import bigfile
import os
import numpy
from multiprocessing import Pool, cpu_count

DefaultHeaderDtype = [
        ('Npart', ('u4', 6)),
        ('Massarr', ('f8', 6)),
        ('Time', ('f8')),
        ('Redshift', ('f8')),
        ('FlagSfr', ('i4')),
        ('FlagFeedback', ('i4')),
        ('Nall', ('u4', 6)),
        ('FlagCooling', ('i4')),
        ('NumFiles', ('i4')),
        ('BoxSize', ('f8')),
        ('Omega0', ('f8')),
        ('OmegaLambda', ('f8')),
        ('HubbleParam', ('f8')),
        ('FlagAge', ('i4')),
        ('FlagMetals', ('i4')),
        ('NallHW', ('u4', 6)),
        ('flag_entr_ics', ('i4')),
    ]

def write_block(block, ff):
    b = block.size * block.dtype.itemsize

    assert b < 2 * 1024 * 1024 * 1024 # avoid overflow!

    b = numpy.array(b, dtype='i4')
    b.tofile(ff)
    block.tofile(ff)
    b.tofile(ff)

def write_gadget_1_ic(filename, header, pos, vel, id):

    with open(filename, 'wb+') as ff:
        write_block(pad256(header), ff)
        write_block(pos, ff)
        write_block(vel, ff)
        write_block(id, ff)

def make_gadget_header(header):
    attrs = header.attrs
    gadget = numpy.zeros((), dtype=DefaultHeaderDtype)
    gadget['Time']  = attrs['Time']
    gadget['Redshift']  = 1.0 / attrs['Time'] - 1
    gadget['Nall'] = numpy.uint32(attrs['TotNumPart'])
    gadget['NallHW'] = numpy.uint32(attrs['TotNumPart'] >> 32)
    gadget['BoxSize'] = attrs['BoxSize'][0] #Mpc/h
    gadget['HubbleParam'] = attrs['HubbleParam']
    gadget['Omega0'] = attrs['Omega0']
    gadget['OmegaLambda'] = attrs['OmegaLambda']
    gadget['Massarr'] = attrs['MassTable']
    print(gadget['BoxSize'])

    return gadget

def pad256(gadget):
    hdtype_padded = numpy.dtype([
                                 ('header', gadget.dtype),
                                 ('padding', ('u1', 256 - gadget.dtype.itemsize))])
    data = numpy.zeros((), dtype=hdtype_padded)
    data['header'] = gadget
    assert data.dtype.itemsize == 256
    return data

def Convert(input, output, Nfile, precision):

    f = bigfile.File(input)
    ds = bigfile.Dataset(f['1/'], ['Position', 'Velocity', 'ID'])

    header = f['Header']
    print("---input file : %s -----", input)
    for key in header.attrs:
        print(key, header.attrs[key])

    gadget_header = make_gadget_header(header)

    dirname = os.path.dirname(os.path.abspath(output))
    if not os.path.exists(dirname):
        print("making dir")
        os.makedirs(dirname)

    exec_convert(output, gadget_header, ds, Nfile, precision)

def exec_convert(basename, baseheader, ds, Nfile, precision):

    print('total number of dm particles', ds.size)
    with Pool(cpu_count()) as pool:
        for i in range(Nfile):
            pool.apply_async(__exec_convert_sep, args=(basename, baseheader, ds, Nfile, precision, i))

def __exec_convert_sep(basename, baseheader, ds, Nfile, precision, i):
    baseheader['NumFiles'] = Nfile
    a = baseheader['Time']

    print('working on file %d/%d' % (i, Nfile))

    start = i * ds.size // Nfile
    end = (i + 1) * ds.size // Nfile

    # read
    data = ds[start:end]
    pos = data['Position'] #Mpc/h
    pos = numpy.array(pos, dtype=precision)

    # convert to gadget 1 units from pecuiliar velocity
    # FIXME: check if this is right.
    vel = data['Velocity'] * a ** -0.5
    vel = numpy.array(vel, dtype=precision)
    id = data['ID']

    filename = '%s.%d' % (basename, i)

    header = baseheader.copy()
    header['Npart'][1] = end - start
    print("header", header)
    write_gadget_1_ic(filename, header, pos, vel, id)
