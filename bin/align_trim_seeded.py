#!/usr/bin/env python3

"""
A wrapper around artic align_trim 1.7.4 that adds in a set seed for consistency
 as we cannot upgrade to the current version at this time with how some of the
 primer schemes are formatted compared to the reference sequence real start and end
"""

import random
from artic.align_trim import main

SEED_VALUE = 42
random.seed(SEED_VALUE)

if __name__ == "__main__":
    main()
