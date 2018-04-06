"""Run all the benchmarks with specific parameters"""
import argparse
import time
from scvi.dataset import load_dataset
from scvi.benchmark import run_benchmarks

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--epochs", type=int, default=250)
    parser.add_argument("-d", "--dataset", type=str, default="cortex")

    # Might be a useful options to combine multiple datasets as an argument
    # parser.add_argument("-l", "--list", help="A list of args", nargs='+', default=[])
    # parser.add_argument("-b", "--bool", help="a bool", action="store_true", default=False)

    args = parser.parse_args()
    gene_dataset = load_dataset(args.dataset)
    start = time.time()
    run_benchmarks(gene_dataset, n_epochs=args.epochs)
    end = time.time()
    print("Total runtime for " + str(args.epochs) + " epochs is: " + str((end - start))
          + " seconds for a mean per epoch runtime of " + str((end - start)/args.epochs) + " seconds.")