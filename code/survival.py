import sys

## Survival data
survival = sys.argv[1]

## Clusters
clust = sys.argv[2]

## Out file: survial_time | event | age | group
outf = sys.argv[3]

## Reading clinical data
patients = {}
f = open(survival,'r')
f.readline()
for line in f:
    pat,time,death,birth = line.strip().split()
    patients[pat] = (time,death,birth)
f.close()

## Reading Cluster data
f = open(clust,'r')
fWrite = open(outf,'w')
print>>fWrite,"Surv_time\tEvent\tAge\tGroup"
for line in f:
    pat,c = line.strip().split()
    print>>fWrite,patients[pat][0],'\t',patients[pat][1],'\t',patients[pat][2],'\t',c
f.close()
fWrite.close()
