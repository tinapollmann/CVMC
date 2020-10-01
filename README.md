This software is used to perform a Monte Carlo simulation of the COVID-19 outbreak dynamics.
It was developed in the context of the ContacTUM collaboration at the Technical
University of Munich.

It requires the root software libraries (https://root.cern) to be installed. 

To build, first edit the makefile, putting in the path to ROOT on your system.

By default, the output consists of a single root file that contains two trees:
1) fPopulationLevelInformation  is ordered by day, and has information on the number of people
exposed, infected, traced, etc .. for each day. If several simulation runs were done,
the trees are all concatenated, ie 'day' is not unique in the tree
2) settings : contains a tree with the values for all the settings

If debug output is enabled, detailed information is printed to std out, a tsv file containing
the same information as the root output file is written, and a .dot output file is
generated, which can be rendered in the graphviz application to create a graphical
view of the infection chain. 


Licensed under MIT Licence https://opensource.org/licenses/MIT
Copyright 2020 ContacTUM
Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without 
restriction, including without limitation the rights to use, copy, modify, merge, publish, 
distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom 
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
