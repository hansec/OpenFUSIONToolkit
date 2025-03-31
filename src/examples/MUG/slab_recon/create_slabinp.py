template = """### Autogenerated slab input file:
## Lx = {0}
## Ly = {1}
## Lz = {2}
vertex 21 xyz 0.0 0.0 0.0 size def * {3}
vertex 22 xyz {0} 0.0 0.0 size def * {3}
vertex 23 xyz 0.0 {1} 0.0 size def * {3}
vertex 24 xyz {0} {1} 0.0 size def * {3}
vertex 25 xyz 0.0 0.0 {2}
vertex 26 xyz {0} 0.0 {2}
vertex 27 xyz 0.0 {1} {2}
vertex 28 xyz {0} {1} {2}

curve 51 vertex 21 22 output yes
curve 52 vertex 23 24 output yes
curve 53 vertex 25 26 output yes
curve 54 vertex 27 28 output yes

curve 61 vertex 21 23 output yes
curve 62 vertex 22 24 output yes
curve 63 vertex 25 27 output yes
curve 64 vertex 26 28 output yes

curve 71 vertex 21 25 output yes
curve 72 vertex 22 26 output yes
curve 73 vertex 23 27 output yes
curve 74 vertex 24 28 output yes

surface 51 curve 61 71 63 73 output yes
surface 52 curve 62 72 64 74 output yes

surface 53 curve 51 71 53 72 output yes
surface 54 curve 52 73 54 74 output yes

surface 55 curve 51 61 52 62 output yes
surface 56 curve 53 63 54 64 output yes

region 1 size def boundary surface -51 52 53 -54 -55 56
"""
#
import sys
nargs = len(sys.argv)
print(nargs)
if nargs != 5:
	print('ERROR: must provide 4 arguments')
	print('Example: ./create_slabinp.py 2.5 3. 3. 1.')
	sys.exit(1)
#
dims = []
for (i,arg) in enumerate(sys.argv):
	if i==0:
		continue
	dims.append(float(arg))
#
print(len(dims))
with open('slab2.inp','w+') as fid:
	fid.write(template.format(dims[0]/2.,dims[1]/2.,dims[2]/2.,dims[3]))
