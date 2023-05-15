import pandas as pd
import os

"""
Add the samples metadata at the end of the humann-summary files.
This file is required for human_barplot script.
"""
def getInputs():
    """
    Get the input file names from the user.

    :return: A tuple of the humann summary file name, the metadata file name, and the output file name prefix.
    """
    paFile = input("Humann summary file (.tsv) name: ")
    metaFile = input("Metadata (in .csv format) filename: ")
    outputFile = input("Output file name prefix: ")
    return paFile, metaFile, outputFile

def checkfile(*args):
    """
    Checks if the given files exist.

    :param args: A list of file names.
    :raises FileNotFoundError: If any of the files do not exist.
    """
    for arg in args:
      if not os.path.exists(arg):
          raise FileNotFoundError(f"{arg} file does not exist")

def main():
    """
    The main function
    """
    # Get the file names from the command line
    pa_file_name, meta_file_name, output_file_name=getInputs()

    # Check if file exists
    checkfile(pa_file_name, meta_file_name)

    # Parse data
    abundance = pd.read_table(pa_file_name, header = 0)
    metadata = pd.read_csv(meta_file_name, header = 0)

    # Work with columns
    column_names = list(abundance.columns)
    col1 = column_names.pop(0) # Poping the first element out from the list
    split_names = [name.split("_") for name in column_names]
    suffix = os.path.commonprefix([name[-1] for name in split_names])

    # Modify the metadata
    metadata.iloc[:,0]=metadata.iloc[:,0]+f"_{suffix}" # Adding string to sample names
    metadata.rename(columns={metadata.columns[0]: col1}, inplace = True) # renaming column
    metadata=metadata.T # Transpose
    metadata.reset_index(inplace=True) # Reset index
    meta = metadata.rename(columns=metadata.iloc[0]).loc[1:] # setting column names

    # Adding samples' metadata to humann output
    df = pd.concat([abundance, meta], join="inner")
    df.to_csv(f'{output_file_name}.tsv', sep="\t", index=False) # Writing output
    print(f"Output file name: {output_file_name}.tsv")


if __name__ == "__main__":
    main()


