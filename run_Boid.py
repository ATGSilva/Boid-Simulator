import sys
from Boid import main

if int(len(sys.argv)) == 2:
    main(int(sys.argv[1]))
else:
    print("Usage: {} <THREADS>".format(sys.argv[0]))

