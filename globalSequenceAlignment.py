import sys
import csv


# Function to read and extract sequences from csv
def csvReader(filename):
    sequences = []

    # Open file and store all sequences in list
    with open(filename, 'r') as sequencesCsv:
        sequencesReader = csv.reader(sequencesCsv)
        next(sequencesReader)  # Skip the first line (header)

        for row in sequencesReader:
            if len(row) >= 2:
                # Join the two sequences with a comma
                sequence = ",".join(row)
                sequences.append(sequence)

    return sequences



# Function to create matrix and call backtracking function to generate result
def globalSequenceAlignment(sequence):
    # Split sequence string into two and define gap penalty
    s1 = sequence.split(',')[0]
    s2 = sequence.split(',')[1]
    gapPenalty = -2

    # Initialize the matrix
    matrix = []
    for j in range(len(s2)+1):
        matrix.append([0]*(len(s1)+1))

    # Initialize first row and column of matrix
    for i in range(len(matrix[0])):
        matrix[0][i] = i*gapPenalty

    for j in range(len(matrix)):
        matrix[j][0] =j*gapPenalty

    # Populate the matrix
    for j in range(len(matrix)):
        for i in range(len(matrix[j])):
            if i!=0 and j!=0:
                # Calculate score using the left and above cells
                gp1 = matrix[j][i-1] + gapPenalty
                gp2 = matrix[j-1][i] + gapPenalty

                # Calculate score for match/mismatch using diagonal cell
                scoringMatrix = matrix[j-1][i-1]
                if s1[i-1]==s2[j-1]:
                    scoringMatrix += 1
                else:
                    scoringMatrix -= 1
                
                # Maximum score among left, above and diagonal is chosen to fill cell
                if gp1>gp2 and gp1>scoringMatrix:
                    matrix[j][i] = gp1
                elif gp2>gp1 and gp2>scoringMatrix:
                    matrix[j][i] = gp2
                else:
                    matrix[j][i] = scoringMatrix
    
    # Print result to console
    print(backtracking(matrix, s1, s2) + " " + str(matrix[-1][-1]))



# Function to perform the backtracking
def backtracking(matrix, s1, s2):
    s1Result = ""
    s2Result = ""
    currPos = [len(matrix)-1, len(matrix[0])-1]

    while currPos[0]>=0 and currPos[1]>=0:
        value = matrix[currPos[0]][currPos[1]]

        # End backtracking when top-left corner of matrix is reached
        if currPos[0]-1<0 and currPos[1]-1<0:
            break

        # Find next position for backtracking, save current values to resulting sequences, and move to next cell
        if value+2 == matrix[currPos[0]-1][currPos[1]]:
            s1Result = "-" + s1Result
            s2Result = s2[currPos[0]-1] + s2Result
            currPos[0] -= 1
        elif value+2 == matrix[currPos[0]][currPos[1]-1]:
            s1Result = s1[currPos[1]-1] + s1Result
            s2Result = "-" + s2Result
            currPos[1] -= 1
        elif value == matrix[currPos[0]-1][currPos[1]-1] + (1 if s1[currPos[1]-1] == s2[currPos[0]-1] else -1):
            s1Result = s1[currPos[1]-1] + s1Result
            s2Result = s2[currPos[0]-1] + s2Result
            currPos = [currPos[0]-1, currPos[1]-1]

    return s1Result + " " + s2Result



# Main function to run everything
def main():
    if len(sys.argv) > 1:
        filename = sys.argv[1]  # Get CSV file from command line argument
        sequences = csvReader(filename)
        
        for sequence in sequences:
            if "," in sequence:
                globalSequenceAlignment(sequence)

# This condition prevents the script from running unless it's executed directly
if __name__ == "__main__":
    main()

