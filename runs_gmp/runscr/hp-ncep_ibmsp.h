#!/bin/sh
#
# @ job_type = parallel
# @ output = out.$(jobid)
# @ error = err.$(jobid)
# @ total_tasks = 12
# @ node = 1
# @ preferences = Feature == "dev"
# @ network.MPI=switch,shared,us
# @ class = dev
# @ wall_clock_limit = 24:00:00
# @ queue
#

MP_SHARED_MEMORY=yes
export MP_SHARED_MEMORY
