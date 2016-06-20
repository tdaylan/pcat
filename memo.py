from memory_profiler import profile
import pcat
from pcat import cnfg

@profile
def main():
    pcat.cnfg.cnfg_ferm_expr_ngal()

if __name__ == '__main__':
    main()
