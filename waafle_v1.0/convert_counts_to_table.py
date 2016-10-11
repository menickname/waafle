#!/usr/bin/python

import argparse
import numpy as np
import re

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def get_args():
    """
    Get arguments passed to script
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="output from count_events table"
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does the following:
    1)
    """
    args = get_args()
    
    dict_matrix = {}
    dict_levelbugavg = {}
    for astrline in open( args.input ):
        aastrline = astrline.strip().split('\t')
        if not re.search( 'header', aastrline[0] ):
            bug1, bug2 = aastrline[0].split('<')
            level = bug1.split('|')[-1].split('__')[0]
            average = np.mean( np.array( [float(x) for x in aastrline[1:]] ) )
            if bug1 == bug2:
                dict_levelbugavg.setdefault( level, [] ).append( [bug1, average] )
            if bug1 in dict_matrix.keys():
                bug2_dict = dict_matrix[bug1]
                bug2_dict[ bug2 ] = average
                dict_matrix[bug1] = bug2_dict
            else:
                dict_matrix[bug1] = {bug2: average} 
    
    for level in dict_levelbugavg.keys():
        fh_matrix = open( level + '_matrix.txt', 'w' )
        sorted_list = sorted( dict_levelbugavg[level], key=lambda bug: bug[1], reverse=True )
        fh_matrix.write( '\t'.join( ['#header'] + [x for x,y in sorted_list] ) + '\n' )
        for bug1, avg1 in sorted_list:
            line = [ bug1 ]
            for bug2, avg2 in sorted_list:
                bug1_dict = dict_matrix[ bug1 ]
                if bug2 not in bug1_dict.keys():
                    value = 0
                else:
                    value = dict_matrix[bug1][bug2]
                line.append( value )
            fh_matrix.write( '\t'.join( str(x) for x in line ) + '\n' )
        fh_matrix.close()
                
    


if __name__ == "__main__":
    main()
