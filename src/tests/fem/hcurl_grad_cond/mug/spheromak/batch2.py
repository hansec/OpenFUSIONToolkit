# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import runsph
import resource
import time
cases = [('mg_bj_np8','bj',8,True),('mg_bjov_np1','bjov',1,True),('mg_cell0_np8','cell0',8,True),
         ('mg_cell1_np32','cell1',32,True),
         ('1lev_bj_np1','bj',1,False),('1lev_bj_np8','bj',8,False),('1lev_bjov_np1','bjov',1,False),
         ('1lev_cell0_np8','cell0',8,False),('1lev_cell1_np32','cell1',32,False)]
for tag, sm, npart, mg in cases:
    t = time.time()
    its, ok = runsph.run(tag, sm, npart, mg)
    print('%-18s ok=%s its=%s wall=%.0f maxrss_GB(max so far)=%.2f' % (tag, ok, its, time.time()-t,
          resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1e6), flush=True)
