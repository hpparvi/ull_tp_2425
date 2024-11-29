
# Timed tests:
These tests were made in the CCA machines, which have 6 cores

For 10 particles, 10k steps:
- no parallel: 0.28 seconds
- parallel: 80.9 seconds

For 100 particles, 10k steps:
- no parallel: 10.1 seconds
- parallel: 197.4 seconds

For 1000 particles, 0.5 seconds:

- without parallelization: 7.6 seconds
- with parallelization: 2.2 seconds

For 1000 particles, 1k steps:
- parallelized: 3.8 seconds
- not parallelized: 14.3 seconds
