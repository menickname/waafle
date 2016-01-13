import os, sys
import waafle_utils as wu

def quote ( string ):
    return "'%s'" % string

for line in open( sys.argv[1] ):
    items = line.strip().split( "\t" )
    items = map( quote, items )
    command = [
        "python",
        "mksyncontig.py",
        "--donor", items[1], 
        "--recipient", items[2], 
        "--donortaxa", items[3], 
        "--reciptaxa", items[4],
		"--taxadiff", items[0],
        "--ngenes", sys.argv[2],
        ]
    os.system( " ".join( command ) )
