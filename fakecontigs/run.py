import os, sys

def quote ( string ):
    return "'%s'" % string

for line in open( sys.argv[1] ):
    items = line.strip().split( "\t" )
    items = map( quote, items )
    command = [
        "python",
        "fakemake.py",
        "--donor", items[1], 
        "--recipient", items[2], 
        "--donortaxa", items[3], 
        "--reciptaxa", items[4],
        "--ngenes", "10",
        ]
    os.system( " ".join( command ) )
