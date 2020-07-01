import os
import strutils

var line=""


let f = open(paramStr(1))

while f.readLine(line):
 echo(line.split()[0])
 echo(line)

close(f)
