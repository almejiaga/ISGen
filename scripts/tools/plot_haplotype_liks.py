import sys, os
import tables
import matplotlib.pyplot as plt
import seaborn

try:
    assert len(sys.argv) == 3
    control_lik_file = os.path.expanduser(sys.argv[1])
    plot_file = os.path.expanduser(sys.argv[2])
except:
    print("Usage: python plot_haplotype_liks.py CONTROL_LIK_FILE PLOT_OUTPUT_FILE")
    sys.exit()

with tables.open_file(control_lik_file, 'r') as f:
    hap_liks = f.root.haplotype_liks[:]

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(sorted(hap_liks))
plt.xlabel('Simulated inheritance path', fontsize=16)
plt.ylabel(r'$log_2$ haplotype likelihood', fontsize=16)
plt.tight_layout()
plt.savefig(plot_file)
