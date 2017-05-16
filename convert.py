import math
import random
import re
import readline
import glob


def convertRowToString(dataRow):
    valuePairStrings = ["{" + str(a) + "," + str(b) + "," + str(c) + "}" for (a, b, c) in dataRow]
    rowString = "{" + ",".join(valuePairStrings) + "}\n"
    return rowString
    
def writeDataToFile(filename, data):
    f = open(filename, 'w')
    for row in data:
        f.write(convertRowToString(row))


def parseFile(filename):
    f = open(filename, 'r')
    data = [parseLine(line) for line in f]
    return data


def parseLine(line):
    inParens = re.search("{(.+)}", line).group(1) #strip the outer parens
    #now we have "[value], {{[value],[value]}, ...}"
    
    rowValValues = re.search("(.+?),.+?{(.+)}", inParens)
    rowValString = rowValValues.group(1)
    rowVal = float(rowValString)
    
    rowValues = rowValValues.group(2) #looks like: "{{[value],[value]}, {[value],[value]}, ...}"
    
    valuePairsParens = re.findall("{.+?}", rowValues) #each pair looks like {[value],[value]}
    
    pairMatch = re.compile("{(.+?),(.+?)}")
    valuePairsString = [pairMatch.match(pair).group(1,2) for pair in valuePairsParens]
    
    
    #no try, the whole parseLine will throw if the numbers don't parse:
    valuePairs = [(float(a), float(b)) for (a,b) in valuePairsString]
    
        
    rowValues = [(rowVal, a, b) for (a, b) in valuePairs]    
    
    return rowValues

def getF(data):
    fe = [f for row in data for (x,y,f) in row]
    return fe

def getColumns(data):
    return len(data[0])

def getRows(data):
    return len(data)

def getRowValues(data):
    return [row[0][0] for row in data]

def getColValues(data):
    return [y for (x,y,z) in data[0]]

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

filename = raw_input("File to open: ")
fileData = parseFile(filename)

fe = getF(fileData)

filenameOut = raw_input("File to save: ")

fileOut = open(filenameOut, 'w')


for f in fe:
    fileOut.write(str(f))
    fileOut.write(",")

