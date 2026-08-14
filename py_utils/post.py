import sys 
from pathlib import Path
from helpers import plot_latest_product

data_dir = Path(sys.argv[1])
figs_dir = str(sys.argv[2])


extent = [40.0, -15.0, 65.0, 30.0, 1.]

products = ["CAL", "DNI", "GHI", "CSR"]

for p in products: 
    plot_latest_product(p, data_dir/p, figs_dir, extent)


