#!/usr/bin/env python3

import sys

"""
This script will perform small tests and compare results with an expected outcome.
Aim is to check if functions in the scripts are producing expected results.
To achive this different functions are imported from extract_topNabundance.py

"""

def compare_dataframe(extracted_df, expected_df):

    
    # 1. Check if number of columns are same:
    print("1. Checking number of columns...")
    if len(extracted_df.columns) == len(expected_df.columns):
        sys.stdout.write(f" \u2714 PASS: Same number of columns = {len(expected_df.columns)}\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(f" \u274c FAIL: Different number of columns (extracted, expected): {len(extracted_df.columns)} != {len(expected_df.columns)}\n")
        sys.stdout.flush()
    
    
    # 2. Check if number of rows are same:
    print("2. Checking number of rows...")
    if len(extracted_df.index) == len(expected_df.index):
        sys.stdout.write(f" \u2714 PASS: Same number of rows = {len(expected_df.index)}\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(f" \u274c FAIL: Different number of rows (extracted, expected): {len(extracted_df.index)} != {len(expected_df.index)}\n")
        sys.stdout.flush()
    
    # 3. Check if the total sum is equal:
    ColSum_extracted = extracted_df.sum(axis = 0) # Column's sum extracted
    ColSum_expected = expected_df.sum(axis = 0) # Column's sum expected

    print("3. Checking total sum...")
    if ColSum_extracted.sum() == ColSum_expected.sum():
        sys.stdout.write(f" \u2714 PASS: Total sum = {ColSum_expected.sum()}\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(f" \u274c FAIL: Total sum (extracted, expected) : {ColSum_extracted.sum()} != {ColSum_expected.sum()} \n")
        sys.stdout.flush()

    print("4. Checking if ColSums are same...")
    # 4. Check if the extracted and expected Colsums dataframe are equal:
    if ColSum_extracted.equals(ColSum_expected):
        sys.stdout.write(f" \u2714 PASS: ColSums are same (extracted, expected)\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(f" \u274c FAIL: ColSums are different (extracted, expected)\n")
        sys.stdout.flush()
    
    print("5. Checking if RowSums are same...")
    RowSum_extracted = extracted_df.sum(axis = 1) # Rows sum
    RowSum_expected = expected_df.sum(axis = 1) # Rows sum

    if RowSum_extracted.equals(RowSum_expected):
        sys.stdout.write(f" \u2714 PASS: RowSums are same (extracted, expected)\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(f" \u274c FAIL: RowSums are different (extracted, expected)\n")
        sys.stdout.flush()

    return