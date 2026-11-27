#!/bin/bash
# mkcase.sh <dir> <dt> <nsteps> <si> <pre: base|pb1> [nlocal] [extra namelist line]
set -e
D=$1; DT=$2; NS=$3; SI=$4; PRE=$5; NL=${6:-440}; EXTRA=$7
SRC=$(dirname "$0")/coax_base_np3
mkdir -p $D
cp $SRC/MUG-coax_railgun-tet.h5 $SRC/plot_railgun.py $D/
sed -e "s/^ nsteps=350 / nsteps=$NS /" -e "s/^ dt=3.e-8 / dt=$DT /" -e "s/^ rst_end=350 / rst_end=$NS /" $SRC/oft_coax.in > $D/oft_coax.in
if [ "$PRE" != "base" ]; then sed -i '' -e 's/^ ittarget=40 / ittarget=100000 /' $D/oft_coax.in; fi
if [ "$SI" != "0" ] || [ -n "$EXTRA" ]; then
  python3 - "$D/oft_coax.in" "$SI" "$EXTRA" <<'PY'
import sys
p, si, extra = sys.argv[1:4]
s = open(p).read()
add = ""
if si != "0":
    add += " xmhd_si=%s\n" % si
if extra:
    add += " %s\n" % extra
s = s.replace(" xmhd_mfnk=T", add + " xmhd_mfnk=T", 1)
open(p, "w").write(s)
PY
fi
if [ "$PRE" == "base" ]; then cp $SRC/oft_coax.xml $D/; else
cat > $D/oft_coax.xml <<EOF
<oft>
  <xmhd>
    <den_bc>1,0,0</den_bc>
    <temp_bc>1,0,1</temp_bc>
    <vel_bc>1,1,2</vel_bc>
    <pb_pre type="1">
      <solver type="gmres"><its>8</its><nrits>8</nrits>
        <pre type="block_jacobi"><nlocal>$NL</nlocal><solver type="lu"></solver></pre>
      </solver>
    </pb_pre>
  </xmhd>
</oft>
EOF
fi
cp $SRC/run.sh $D/; cp $SRC/case.json $D/
