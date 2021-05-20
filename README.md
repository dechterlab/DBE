# DBE

Deep Bucket Elimination based on code by Kalev.

## Instructions

Do the cmake for the first time to have the cmakefiles.
```

cd BESampling
mkdir build
cd build
cmake .. -B . -DCMAKE_PREFIX_PATH=~/Downloads/libtorch  //path to libtorch

```

```

make

./BESampling -fUAI ../3testproblems/grid20x20.f5.uai -fVO ../grid20x20.f5.uai.ord.elim -iB 9999 -v2sample 369 -nsamples 10 -h_dim 100 --network net

./BESampling -fUAI ../3testproblems/pedigree19.uai -fVO ../3testproblems/pedigree19.uai.ord.elim -iB 9999 -v2sample 369 -nsamples 10 -h_dim 100  --network masked_net
```

### Sample Result
```
STATS : nBuckets=387 MaxNumVarsInBucket=20 nPartitionedBuckets=0
wMBE done; result=-29.6304 runtime=35sec tablesmemory=534828304 bytes
largest bucket v=369
will sample v=369 nsamples=10 bSize=20 vars nMBs=1
output fn size=1769472 entries
fn signature : 74 75 80 97 99 215 217 224 226 247 271 297 319 362 363 364 98 218 361
sample : 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 2 2 0 0 = -inf logscale (idx=842540)
sample : 1 0 0 1 1 0 1 0 1 0 0 1 0 0 1 0 2 1 1 = -inf logscale (idx=1068383)
sample : 1 0 0 0 1 0 0 0 1 0 1 0 0 0 0 2 1 1 0 = -inf logscale (idx=944382)
sample : 1 1 0 1 0 0 1 1 0 1 1 0 1 1 0 1 1 1 0 = -29.7527 logscale (idx=1461330)
sample : 0 0 0 1 1 1 0 1 0 0 0 1 0 1 0 0 2 1 0 = -inf logscale (idx=200962)
sample : 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 2 0 0 1 = -inf logscale (idx=897721)
sample : 1 0 0 0 0 0 0 1 0 1 0 0 1 2 0 1 2 1 1 = -inf logscale (idx=893759)
sample : 1 1 1 1 0 0 1 1 0 1 1 0 0 1 0 2 0 0 1 = -inf logscale (idx=1682305)
sample : 0 1 0 1 0 1 0 0 1 1 1 0 0 2 0 2 2 1 1 = -22.6421 logscale (idx=586835)
sample : 1 1 0 1 1 1 0 1 0 0 0 1 0 2 1 2 2 1 1 = -25.5474 logscale (idx=1528199)
```
