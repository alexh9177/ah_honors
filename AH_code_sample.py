

# Alex Hirano
# Biofrontiers Institute 

# Transcriptome analysis pipeline

# Code setup and imports
import re
from collections import Counter
import csv

def SalmonInput(salmonFile):
    '''

    TO-DO -> find salmon quant header data

    Returns a nested list of Salmon quantified-isoforms.

        Parameters:
            salmonFile (str): Absolute path to the 'quant.sf' Salmon quantification file.
                EX: "/path/to/quant.sf"

        Returns:
            salmonList (nested list): List of isoforms by stringtie header w/ 
            Salmon quantification data.
                Each nested list element contains:
                - 
                - 
                - 
                - 
    
    '''
    salmonList = []
    with open(salmonFile) as salmonFile:
        # Removes header
        next(salmonFile) 
        for line in salmonFile:
            salmonList.append(line.split()) 
    return salmonList

def BlastInput(blastFile):
    '''
    Returns a filtered, nested list of BLASTn isoform matches.

        Parameters:
        -----------
            blastFile (str): String of absolute path to the BLASTn output file.
                EX: "/path/to/blastn_output.txt"
    
        Returns:
        --------
            blastList (nested list): List of all isoforms w/BLASTn data.
                Each nested list element contains:
                - Stringtie name (from hybrid transcriptome)
                - Closest human Ensembl transcript match
                - Closest human gene homolog name
                - Gene type
    '''
    blastList = []
    with open(blastFile) as blastFile:
        for line in blastFile:
            # Uses regex to remove blastn pseudogene hits
            # Catches all variants of "pseudogene, processed_pseudogene, unprocessed_pseudogene"
            pseudoSort = re.search("pseudogene",line)
            if not pseudoSort: 
                # Command line blastn output has "|" char between some variables
                # Replaces line char with space before split
                tempLine = line.replace("|", " ")
                tempList = tempLine.split()
                # Selects stringtie name [0], ensembl transcript [1], gene name [6]
                # and gene type [8], appends to blastList.
                blastList.append(
                    [tempList[0],tempList[1],
                    tempList[6],tempList[8]]
                    )
        return blastList


def ReadEncode(encodeFile):
    '''
    Returns a list of Encode transcripts from a downloaded text file.

        !! IMPORTANT !!
        --------------- 
            - Text file should not contain any headers
            - Each transcript name should be on its own line
            - Otherwise strange things may happen

        Parameters:
        -----------
            encodeFile (str): String of absolute path to the BLASTn output file.
                EX: "/path/to/human_immune_transcripts.txt"
    
        Returns:
        --------
            encodeSortList (list): List of all encode transcripts in the text file.
    '''
    encodeSortList = []
    with open() as sortFile:
        for line in sortFile:
            # Sometimes newline characters cause issues
            encodeSortList.append(line.strip('\n')) 

    return encodeSortList


def AnnotationMerge(salmonList, blastList):
    '''
    Returns a nested list of quantified annotations.

    The functions SalmonInput and BlastInput MUST be used to clean the respective lists beforehand.

        Parameters:
        -----------
            salmonList (nested list): List of isoforms by stringtie header w/ Salmon quantification data.
            blastList (nested list): List of all isoforms w/BLASTn data.
            encodeSortList (list): List of all encode transcripts in the encode text file.

            isoformFile (str): String of absolute path to file location.
                File contains 100 human immune isoforms from ENCODE database.
                EX: "path/to/human_immune_genes.txt"
    
        Returns:
        --------
            annotationList (nested list): Merged Salmon and BLASTn data for each transcript.
    '''
    annotationList = []
    for i in blastList:
        for j in salmonList:
            # Checks for matching stringtie headers
            if i[0] == j[0]:
                # Removes stringtie header from salmon isoform, preventing stringtie duplication
                j.pop(0)  
                # Concatenates the STRG header-matched strings
                annotation = i + j 
                annotationList.append(annotation)

    print(
        "{0} annotations were generated.".format(len(annotationList))
        )

    return annotationList


def AnnotationSort(annotationList, encodeSortList):
    '''
    Returns a filtered, nested list of quantified annotations.

    The functions AnnotationMerge and ReadEncode MUST be used first to generate the annotation and encode
    lists before sorting.

        Parameters:
        -----------
            annotationList (nested list): Merged Salmon and BLASTn data for each transcript.
            encodeSortList (list): List of all encode transcripts in the encode text file.

            isoformFile (str): String of absolute path to file location.
                The test file contains 100 human immune isoforms from ENCODE database.
                EX: "path/to/human_immune_genes.txt"
    
        Returns:
        --------
            sortedAnnotationList (nested list): Merged Salmon and BLASTn data sorted by encodeSortList.         
    '''
    sortedAnnotationList = []    
    # Removes the period from the blastn encode transcript data to search against 
    # list of downloaded transcripts from encode database
    for annotation in annotationList:
        # Splits by period, sets annotation[1] (encode transcript name) to everything before period
        a = annotation[1].split(".")
        b = a[0]
        annotation[1] = b

    # Checks if encode transcript matches annotation, appends to sortedAnnotationList if true
    sortedAnnotationList = []
    for annotation in annotationList:
        for transcript in encodeSortList:
            if annotation[1] == transcript:
                sortedAnnotationList.append(annotation)

    print(
        "{0} sorted annotations were generated.".format(len(sortedAnnotationList))
        )

    return sortedAnnotationList

def IsoformCounter(aList, sortFlag = None):
    '''
    Returns a dictionary of isoforms isoformDict, where the key-value pair is isoform:frequency.

        Parameters:
        -----------
            aList (nested list): Generated from AnnotationSort, will take sorted/unsorted inputs.
            sortFlag (boolean): Set to True to print isoformDict from high isoform count -> low.
    
        Returns:
        --------
            isoformDict (dict): Frequency count of isoforms.
    '''
    nameList = []
    for isoform in aList:
        # Pulls gene name from all isoforms to make list, then counts
        nameList.append(isoform[2])
    
    # Uses counter to generate dictionary from nameList
    isoformDict = Counter(nameList)

    # Sort function enabled 
    if sortFlag == True:
        # Uses most_common and dictionary length to sort dictionary
        isoformDict = isoformDict.most_common(len(isoformDict))

    return isoformDict


def NestedLstCSVWriter(nestedList, csvHeaderList, fileName = None):
    '''
    Writes a nested list to a csv file, where each csv line is a list element.

        Parameters:
        -----------
            nestedList (nested list): Any nested list of strings.
            headerList (list): A list of strings to use as the CSV header.
            fileName (str): Any string to name the csv output.
                Output will be formatted as "fileName_annotation.csv"
                If no fileName set, default format is "output_annotation.csv"
    
        Returns:
        --------
            None
    '''
    # Error checker for CSV header
    if all(isinstance(i, str) for i in csvHeaderList) == False:
        print("ERROR: csvHeaderList must be a list of strings.")

    # Error checker for fileName
    # Handles conditions where no name entered
    if fileName == None:
        fileName = "output_annotation.csv"
    elif type(fileName) != str:
        print("ERROR: Please enter a string for the file name.")
        return 
    else:
        # Adds an underscore for file naming
        fileName = fileName + "_annotation.csv"
    # Begin writer 
    outFile = open(fileName, "w")
    writer = csv.writer(outFile)
    # Write the header first
    writer.writerow(csvHeaderList)

    # Iterate through nested list and write each line
    for list in nestedList:
        writer.writerow(list)
    return 


