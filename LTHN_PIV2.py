# -*- coding: utf-8 -*-
#"""
#
# (               )     )   |
# )\ )  *   )  ( /(  ( /(   |
#(()/(` )  /(  )\()) )\())  |
# /(_))( )(_))((_)\ ((_)\   | LTHN PIV: Is an opensource PIV toolbox
#(_)) (_(_())  _((_) _((_)  |
#| |  |_   _| | || || \| |  |
#| |__  | |   | __ || .` |  | Website: https://github.com/esterfiorillo/PIV-method-implementation-code
#|(___|(|_|   |_||_||_|\_|  |
#                           | 
# )\ ) )\ )                 | CDTN - Centro de Desenvolvimento da Tecnologia Nuclear
#(()/((()/( (   (           | LTHN - Laboratório de Termo-Hidráulica e Neutrônica
# /(_))/(_)))\  )\          | Belo Horizonte, MG, Brasil
#(_)) (_)) ((_)((_)         |
#| _ \|_ _|\ \ / /          | @authors: esterfiorillo, acampagnole 
#|  _/ | |  \ V /           |
#|_|  |___|  \_/            |
#
#
#"""
# License
#     This file is part of LTHN PIV.
#
#     LTHN PIV is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     LTHN PIV is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with LTHN PIV.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     ./LTHN2.py

"""
Only this file needs to be executed to generate the output by using input txt files and not the interface.
"""
"""
   Load basic modules
"""
import sys

"""
    try to load processing functions
"""
try:
    from processing.read_input_files import processing_input_file
except Exception as E:
    print("Missing critical processing module. Please install module:")
    print(E)
    sys.exit()
    

def main():
    """
    Reads the string that has the address of the input file and call the function that executes the process. 
    """
    f_address = 'Input_files/file1'
    
    processing_input_file(f_address)

if __name__ == "__main__":
    main()