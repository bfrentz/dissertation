#!/usr/bin/env python3

import os, sys

f = open('./chapter1.tex','r')
filedata = f.read()
f.close()

newdata = filedata.replace(r"\\cite{","\cite{")

f = open('./chapter1.tex','w')
f.write(newdata)
f.close()
