#!/usr/bin/env python3
import urllib.request
import sys
import os

import spaTyper.enricher
import spaTyper.utils

####################################################
def getSpaTypes(reps, orders):
    """
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    
    """
    seqDict = {}
    letDict = {'58': 'B4', '30': 'O2', '54': 'H3', '42': 'M2', '48': 'V2', '45': 'A3', '43': 'X2',
               '60': 'S2', '61': 'W3', '62': 'U3', '57': 'S', '64': 'X3', '49': 'Y2', '66': 'F4',
               '90': 'I', '68': 'E4', '69': 'C4', '80': 'K4', '52': 'R3', '53': 'G3', '02': 'A',
               '03': 'D2', '26': 'T', '01': 'XX', '06': 'G2', '07': 'U', '04': 'Z', '05': 'C',
               '46': 'Y3', '47': 'Z3', '08': 'X', '09': 'A2', '28': 'R', '29': 'F2', '41': 'U2',
               '14': 'I2', '59': 'T3', '78': 'J4', '51': 'P2', '24': 'Q', '56': 'J2', '25': 'O',
               '39': 'E3', '65': 'S3', '76': 'K3', '75': 'I4', '38': 'F3', '73': 'G4', '72': 'P3',
               '71': 'Q3', '70': 'D4', '20': 'D', '74': 'H4', '21': 'F', '11': 'Y', '10': 'C2',
               '13': 'E', '12': 'G', '15': 'W', '22': 'L', '17': 'M', '16': 'K', '19': 'H', '18': 'H2',
               '31': 'N', '23': 'J', '37': 'D3', '36': 'W2', '35': 'C3', '34': 'B', '33': 'P', '55': 'A4',
               '63': 'V3', '32': 'E2', '44': 'Z2', '50': 'T2'}
    typeDict = {}
    seqLengths = set()

    reps_dict = spaTyper.utils.fasta_dict(reps)
    for i in reps_dict:
        seq = reps_dict[i]
        num = i[1:]
        seqDict[seq.upper()] = num
        seqLengths.add(len(seq))
    with open(orders) as f:
        for line in f:
            st, pattern = line.rstrip().split(',')
            typeDict[pattern] = st
    return seqDict, letDict, typeDict, seqLengths

####################################################
def findPattern(infile, seqDict, letDict, typeDict, seqLengths, enrich):
    """
    Finds the Spa type given the repeat order
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    """
    qDict = spaTyper.utils.fasta_dict(infile)
    if enrich:
        seq_list = spaTyper.enricher.enrichSeq.check_primers()
    else:
        seq_list = []
        for i in qDict:
            seq_list.append(qDict[i])
    rep_list = []
    for i in seq_list:
        index = 0
        adjacent = False
        rep_order = []
        while index <= len(i):
            gotit = False
            for j in seqLengths:
                if i[index:index+j] in seqDict:
                    if adjacent or rep_order == []:
                        rep_order.append(seqDict[i[index:index+j]])
                    else:
                        rep_order.append('xx')
                        rep_order.append(seqDict[i[index:index+j]])
                    index += j
                    gotit = True
                    adjacent = True
                    break
            if not gotit:
                index += 1
                adjacent = False
        rep_list.append(rep_order)
    out_list = []
    for i in rep_list:
        let_out = ''
        for j in i:
            if j in letDict:
                let_out += letDict[j] + '-'
            else:
                let_out += 'xx-'
        let_out = let_out[:-1]
        if '-'.join(i) in typeDict:
            type_out = typeDict['-'.join(i)]
        else:
            type_out = '-'.join(i)
        out_list.append(let_out)
        out_list.append(type_out)
    return out_list
