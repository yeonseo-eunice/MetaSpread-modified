import numpy
import random as rand
import pytest

# We need to set a seed for consistency
# This seed specifically puts no more that 2 cancercells per gridcell,
# allowing to test proliferation
@pytest.fixture
def random():
    rand.seed(1)