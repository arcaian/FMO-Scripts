# Import glob to check for all files in the directory, os to create the directory, and shutil to move files to the directory
import glob
import os
import shutil


# Declare a function to read the input file's fragment atom listing
def extractFragmentBounds(inputFile, startSearch, endSearch):
    inputFile = open(inputFile, 'r')
    readMode = False
    currentFragmentAtoms = []
    fragmentLists = []  # We declare the input file, and a series of empty variables for later
    for line in inputFile:
        if endSearch in line:  # For each line in the input file, stop splitting lines and reading if the endSearch tag is present
            readMode = False
        if readMode is True:
            # If readMode is turned on, split the line at blank space, and take out any purely blank space
            splitLine = [x.strip() for x in line.split()]
            currentAtom = 0
            for atom in splitLine:
                previousAtom = currentAtom
                # For each atom in the split line, keep the previous atom stored, and keep an integer version of the current atom
                currentAtom = int(atom)
                if currentAtom > 0:  # If the atom is over 0, then just append it to the list
                    currentFragmentAtoms.append(currentAtom)
                elif currentAtom < 0:  # But if the atom is under 0, that's FMO terminology for all atoms from the previous one to this one are included; these are all added
                    for i in range((previousAtom+1), (abs(currentAtom)+1)):
                        currentFragmentAtoms.append(i)
                elif currentAtom == 0:  # A 0 marks the end of a fragment; this means we append the list of atoms in the fragment to the total list, and reset the counter
                    fragmentLists.append(currentFragmentAtoms)
                    currentFragmentAtoms = []
        # When the string corresponding to the start of the coordinate section is found, enable readMode (at the end, so the next line is the first split line)
        if startSearch in line:
            readMode = True
    return(fragmentLists)  # When the whole set of loops is done, return the list of all fragments


# A very simple function - it finds the fragment an atom is associated with
def findAtomsFragment(fragmentList, atomNumber):
    fragmentCounter = 0
    for fragment in fragmentList:
        for atom in fragment:  # All atoms in each fragment are looped through, checking for a match with the atom of interest; when found, the fragment number is returned
            if atomNumber == atom:
                return(fragmentCounter)
        fragmentCounter += 1


def findBondedAtoms(fragmentList, inputFile):
    inputFile = open(inputFile, 'r')
    readMode = False
    bdaSetList = []
    for line in inputFile:
        if '$END' in line:  # For each line in the input file, stop splitting lines and reading if the endSearch tag is present
            readMode = False
        if readMode is True:
            # If readMode is turned on, split the line at blank space, and take out any purely blank space
            splitLine = [x.strip() for x in line.split()]
            BDA = abs(int(splitLine[0]))
            BAA = int(splitLine[1])
            setToAdd = [BDA, BAA]
            # We don't want the '-' that's present to indicate its on the previous fragment, so take the absolute value, and then add the BDA/BAA to a list of lists
            bdaSetList.append(setToAdd)
        # When the string corresponding to the start of the coordinate section is found, enable readMode (at the end, so the next line is the first split line)
        if '$FMOBND' in line:
            readMode = True
    fragmentCounter = 0
    newFragmentList = fragmentList
    for fragment in newFragmentList:
        for bdaSet in bdaSetList:  # Loop through each of the sets of BDA/BAA
            BDA = 'F' + str(bdaSet[0])
            # We then define the BDA/BAA as strings with 'F' at the start, for fragment - so they're differentiable later for conversion to hydrogens
            BAA = 'F' + str(bdaSet[1])
            for atom in fragment:
                # For each atom, we check if it matches the BAA; if it does, we append the BDA to this fragment (as it's normally part of the previous fragment)
                if bdaSet[1] == atom:
                    fragment.append(BDA)
                    # Wew then find the BDA's fragment, and append the BAA to it, as it's not normally part of that fragment
                    BDAFrag = findAtomsFragment(fragmentList, bdaSet[0])
                    newFragmentList[BDAFrag].append(BAA)
        # Once done with each of the BDA sets, append the fragment to the growing list of fragments and increment a counter for the total amount of fragments dealt with
        fragmentCounter += 1
    return(newFragmentList)


# Now we need a function to take those atom numbers and turn them into individual fragment .xyz files
def fetchAtomCoord(inputFile, atomLists, startSearch, endSearch):
    inputFile = open(inputFile, 'r')
    readMode = False
    atomDictionary = {}
    residueSet = []  # Again, the input file is opened to read, and variables are initialized
    for line in inputFile:
        if endSearch in line:  # The same logic as before; stop reading when you find the endSearch string
            readMode = False
        if readMode is True:
            # While reading, split the line and strip whitespace out
            splitLine = [x.strip() for x in line.split()]
            # The current atom are all numbers after the first character in the atom number, which is the element of the atom
            currentAtom = int(splitLine[0][1:])
            # Define the atomic information as the element, plus all coordinate information
            atomInformation = [splitLine[0][0], splitLine[2:]]
            # Append to a dictionary the atomic information, with the atom number as the registry item
            atomDictionary[currentAtom] = atomInformation
        if startSearch in line:
            readMode = True
    for fragment in atomLists:  # With the above loop completed, take the list of fragment atom IDs from the extractFragmentBounds function and loop through it
        fragmentCoords = {}
        for atom in fragment:  # For each atom in the fragment, match the atom number with its atomic information
            # If it's a string, then it's one of the BDA/BAA atoms, so convert it back to an int, and replace the element with hydrogen
            if type(atom) == str and str(atom)[0] == 'F':
                atom = int(atom[1:])
                getValue = atomDictionary.get(atom)
                # To avoid issues with modifying the original dictionary, we make a separate variable with the contents, delete the existing element type, and add in 'H'
                toAddSet = getValue[1:]
                toAddSet.insert(0, 'H')
                # This could potentially be made variable, but hydrogen is the simplest - no need to figure out other capping atoms
            else:  # If it's not a string, it can be treated normally
                toAddSet = atomDictionary[atom]
            fragmentCoords[atom] = toAddSet
        # And append this to a dictionary, giving a list of dictionaries, each with atom number and atomic information correlated
        residueSet.append(fragmentCoords)
    return(residueSet)


# With the residues having their atoms + coordinates associated with them, we now need to write the .xyz files
def makeIndividualFiles(coordDict, outputPrefix):
    fragmentCounter = 0
    # Set the output folder as the current path plus the outfix prefix, and make that folder if it doesn't exist
    outputFolder = os.getcwd() + '/' + outputPrefix + '_xyz'
    if(os.path.exists(outputFolder) is False):
        os.mkdir(outputFolder)
    for fragment in coordDict:  # For each fragment in the coordinate dictionary, write the first two lines - the number of atoms, and a comment
        fragmentCounter += 1
        outputName = outputPrefix + str(fragmentCounter) + '.xyz'
        outputFile = open(outputName, 'w')
        numberOfAtoms = str(len(fragment)) + '\n'
        outputFile.write(numberOfAtoms)
        commentMessage = 'An xyz file generated from fragment ' + \
            str(fragmentCounter) + ' of an FMO input file for GAMESS-US\n'
        outputFile.write(commentMessage)
        for coord in fragment:  # For each coordinate set in the fragment, separate out the x/y/z coordinates and make sure they've got appropriate spacing for a neat file
            currentCoord = fragment[coord][1]
            currentAtom = str(fragment[coord][0])
            xCoord = currentCoord[0]
            while len(xCoord) < 10:
                xCoord = ' ' + xCoord
            yCoord = currentCoord[1]
            while len(yCoord) < 10:
                yCoord = ' ' + yCoord
            zCoord = currentCoord[2]
            while len(zCoord) < 10:
                zCoord = ' ' + zCoord
            # And then set a string of coordinates with the element name at the start, and write that to file
            stringCoord = currentAtom + ' ' + xCoord + yCoord + zCoord + '\n'
            outputFile.write(stringCoord)
        currentFilePath = os.getcwd() + '/' + outputName
        # With the .xyz file written, move it to the output folder - much neater than a big set of files in the working directory
        shutil.move(currentFilePath, (outputFolder + '/' + outputName))
    print('XYZ files have been generated from fragments in the following folder: ' + outputPrefix + '_xyz')


# This is a relatively simple function that takes in a list of fragments, and outputs a valid Psi4 input file for the coordinates
# Currently there's the big issue of these inputs not being capped - and so are unpaired electrons; I'm not certain on how to fix this one
def makePsi4Input(fragmentList, fragNums, inputFile, generatedFileName, basisSet, energyType, memoryRequired):
    fragFileString = ''
    fragComString = ''
    for num in fragNums:  # Simply splitting the given list of fragment numbers into either _ or , deliniated values for use in naming files/adding comments
        fragFileString = fragFileString + '_' + str(num)
        fragComString = fragComString + ', ' + str(num)
    writeFileName = generatedFileName + fragFileString + '.dat'
    psi4File = open(writeFileName, 'w')
    psi4File.write('# This is a Psi4 input file generated from the fragment(s)' +
                   fragComString + ' of the FMO input file ' + inputFile + '\n')
    psi4File.write('# This is an ' + energyType[1:-1] + '/' +
                   basisSet + ' single point energy calculation')
    # A set of two comments for the start of the input file, explaining the source of the geometry and the calculation to be run
    # Specifying the amount of RAM required for the job
    psi4File.write('\nmemory ' + str(memoryRequired) + ' gb \n')
    moleculeString = 'molecule {\n'  # This is the start of the xyz section for the molecule
    for fragment in fragmentList:
        for atom in fragment:  # For each atom in each fragment, we have to do some processing
            elementName = fragment[atom][0]
            xyzList = fragment[atom][1]
            fragment[atom][1][0] = format(float(fragment[atom][1][0]), '.4f')
            fragment[atom][1][1] = format(float(fragment[atom][1][1]), '.4f')
            fragment[atom][1][2] = format(float(fragment[atom][1][2]), '.4f')
            xValue, yValue, zValue = xyzList[0], xyzList[1], xyzList[2]
            # Firstly, we round all the co-ordinates to 4 decimal places and assign them variables
            while len(xValue) < 10:
                xValue = ' ' + xValue
            while len(yValue) < 10:
                yValue = ' ' + yValue
            while len(zValue) < 10:
                zValue = ' ' + zValue
            # Then we make each co-ordinate 10 characters long, filling in blank space; this ensures consistent formatting for positive and negative values so long as they are under 1000
            stringToAdd = '  ' + elementName + ' ' + xValue + ' ' + yValue + ' ' + zValue + '\n'
            moleculeString = moleculeString + stringToAdd
            # Finally, we add all the disparate bits together to make a line for the molecule in Psi4, and add that to a string that represents all the coordinates
    moleculeString = moleculeString + '}\n'
    psi4File.write(moleculeString)
    psi4File.write('\nset basis ' + basisSet)
    psi4File.write("\nenergy(" + energyType + ")")
    # And finally, we write the large molecule string, as well as strings defining the basis set and energy type

# Then we just do a very small function - makes it neater and easier to manage. It calls the above 3 functions, going from input file to written .xyz files.


def makeFragments(inputFile, outputName):
    fragmentList = extractFragmentBounds(inputFile, 'INDAT(1)', 'RESPAP')
    bondedList = findBondedAtoms(fragmentList, inputFile)
    coordList = fetchAtomCoord(inputFile, bondedList, '$FMOXYZ', '$END')
    makeIndividualFiles(coordList, outputName)
    # 3 easy calls to functions that make a list of fragments, fetch their coordinates, and make the .xyz file
    print('The following questions relate to the following GAMESS-US input file:', inputFile)
    fragmentString = input(
        "For which fragments (if any) would you like to create Psi4 input files? Please separate by commas (None): ") or 'None'
    if fragmentString != 'None':
        basisSet = input(
            "What is your preferred basis set, in Psi4 formatting (cc-pVTZ): ") or 'cc-pVTZ'
        energyType = input(
            "What energetic method is to be used, in Psi4 formatting ('mp2'): ") or "'mp2'"
        memoryRequired = input(
            "How many GB of memory will be required for this job (8): ") or "8"
        # We then define some important factors; the defaults are listed in brackets
        fragmentsOfInterest = []
        for fragment in fragmentString.split(','):
            fragmentsOfInterest.append(int(fragment))
            totalFragments = []
        for fragment in fragmentsOfInterest:
            fragmentCoord = coordList[fragment - 1]
            totalFragments.append(fragmentCoord)
            makePsi4Input([fragmentCoord], [fragment], inputFile,
                          outputName, basisSet, energyType, memoryRequired)
            # For each fragment in our list of fragments of interest, we call the makePsi4Input function on the individual fragments, but also collect the coordinate info in a larger list
        makePsi4Input(totalFragments, fragmentsOfInterest, inputFile,
                      outputName, basisSet, energyType, memoryRequired)
        # That larger list is then used to call the makePsi4Input function on the entire set of fragments


# Make a list of all .inp files in the working directory, and for each of them, run the makeFragments function with an output name that is the same as the input file, without the .inp extension name
inputList = glob.glob('*inp')
for file in inputList:
    outputName = file[:-4]
    makeFragments(file, outputName)
