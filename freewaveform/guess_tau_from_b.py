import numpy as np 

def main(target_b, Gmax):
    magic_number = 702129.6736002087
    tau = ((float(target_b)*1e9) / (float(Gmax)**2*magic_number))**(1/3.)
    print(tau)
    return tau

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
