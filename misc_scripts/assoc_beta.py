#!/usr/bin/env python

import os, sys, re
import argparse, textwrap
import pandas as pd
import numpy as np

from math import log
#from scipy.stats import norm

def isLogisticResult(_lr):


    df_NULL = pd.DataFrame([])


    if os.path.exists(_lr):

        if _lr.endswith('.assoc.logistic'):

            # PLINK Logistic Regression Result.
            df_RETURN = pd.read_csv(_lr, sep='\s+', header=0, usecols=['SNP', 'OR', 'SE', 'STAT', 'P'])
            df_RETURN.dropna(inplace=True)

            df_RETURN['BETA'] = df_RETURN['OR'].map(lambda x : log(x))

            return df_RETURN

    else:
        print("Given Logistic Regression result file('{}') doesn't exist.\n"
                                            "Please check it again.".format(_lr))
        return df_NULL


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description = 'Add BETA column to assoc.logistic file.')

    argparser.add_argument("--logistic-result", "-lr", help="\nLogistic Regression Result file (ex. *.assoc.logistic).\n\n", required=True)
    argparser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)

    args = argparser.parse_args()
    print(args)

    isLogisticResult(args.logistic_result).to_csv(args.out, sep='\t', header=True, index=False)



