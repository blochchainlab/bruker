import numpy as np 

def main(target_b, tau):
    magic_number = 702129.6736002087
    G = (((float(target_b)*1e9)/(float(tau))**3) / magic_number)**0.5
    print(G)
    return G

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
