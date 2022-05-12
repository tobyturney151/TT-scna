# Import Modules
from ttscna import ttscna
import pandas as pd
from os.path import exists

# Generate some data
ttscna.ttscna('Example_Data/2021_10_22_0029.atf', 95.0, 260.0)

# Read the data you generated
results = pd.read_csv('Example_Data/2021_10_22_0029_results.csv', index_col = 0)

# Was the .csv made?
def test_csv_exists():
    assert exists('Example_Data/2021_10_22_0029_results.csv') == True

# Was the .png made?
def test_png_exists():
    assert exists('Example_Data/2021_10_22_0029_results.png') == True

# Is the Unitary Conductance what we'd expect?
def test_uconn():
    assert results['Unitary Conductance (fS)'][0] > 400.0 and results['Unitary Conductance (fS)'][0] < 500.0

# Is the Mean Dwell Time what we'd expect?
def test_dwell():
    assert results['Dwell Time (ms)'][0] > 5000.0 and results['Dwell Time (ms)'][0] < 6000.0

# Is the Mean Bulk Current what we'd expect?
def test_subcurr():
    assert results['Mean Bulk Current (pA)'][0] > -150.0 and results['Mean Bulk Current (pA)'][0] < -50.0
