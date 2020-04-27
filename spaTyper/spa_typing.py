#!/usr/bin/env python3
import urllib.request
import sys
import os

import spaTyper.enricher
import spaTyper.utils

####################################################
def getSpaTypes(reps, orders, debug):
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
def findPattern(qDict, seqDict, letDict, typeDict, seqLengths, enrich, debug):
    """
    Finds the Spa type given the repeat order
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    """
    
    ### create list of sequences: either enrich or all sequences
    if enrich:
        if debug:
            print ("## Debug: enrich sequences with primer seqs")
        
        rep_list = []
        seq_list = spaTyper.enricher.enrichSeq.check_primers(qDict)
        for i in seq_list:
            pattern = findPattern_sequence(i, seqDict, seqLengths, debug)
            if pattern:
                type_return = findPattern_type(pattern, letDict, typeDict)
                rep_list.append(type_return)
        
        return (seq_list)
        
    else:
        if debug:
            print ("## Debug: use all sequences")
        
        dict_repeats = {}
        for i in qDict.keys():
            pattern = findPattern_sequence(qDict[i], seqDict, seqLengths, debug)
            if pattern:
                type_return = findPattern_type(pattern, letDict, typeDict)
                dict_repeats[i] = type_return
                
        return (dict_repeats)

####################################################
def findPattern_sequence(seq, seqDict, seqLengths, debug):
    """
    Identify the pattern of repeats per sequence
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    """
    
    index = 0
    adjacent = False
    rep_order = []
    while index <= len(seq):
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

    ## debugging nessages
    if debug:
        print ('## Debug: rep_order:')
        print ("rep_order: ", rep_order)

    ## if it is not empty
    if rep_order:        
        return (rep_order)
    else:
        return()

####################################################
def findPattern_type(pattern, letDict, typeDict):
    """
    Identifies the SPA repeat type
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    
    """
   ## 
    let_out = ''
    for j in pattern:
        if j in letDict:
            let_out += letDict[j] + '-'
        else:
            let_out += 'xx-'
    let_out = let_out[:-1]
    if '-'.join(pattern) in typeDict:
        type_out = typeDict['-'.join(pattern)]
    else:
        type_out = '-'.join(pattern)
        
     ## debugging nessages
    if debug:
        print ('## Debug: rep_list: (let_out, type_out)')
        print ('rep_list pattern:', pattern)
        print ('let_out', let_out)
        print ('type_out', type_out)
            
    string_return = let_out + '::' + type_out
    return (string_return)

####################################################
def OLD_findPattern(infile, seqDict, letDict, typeDict, seqLengths, enrich, debug):
    """
    Finds the Spa type given the repeat order
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    """
    qDict = spaTyper.utils.fasta_dict(infile)

    ### create list of sequences: either enrich or all sequences
    if enrich:
        if debug:
            print ("## Debug: enrich sequences with primer seqs")

        seq_list = spaTyper.enricher.enrichSeq.check_primers()

    else:
        if debug:
            print ("## Debug: use all sequences")
 
        seq_list = []
        for i in qDict:
            seq_list.append(qDict[i])
    
    ## create repeat identification
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

        ## if it is not empty
        if rep_order:        
            ## debugging nessages
            if debug:
                print ('## Debug: seq_list rep_order:')
                #print ("Seq: i", i)
                print ("rep_order: ", rep_order)

            rep_list.append(rep_order)

    ## return if empty
    if not rep_list:
        return('NA')
    
        ## debugging nessages
    if debug:
        print ('## Debug: rep_list:')
        print (rep_list)

    ## Identify the repeat type
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
            
         ## debugging nessages
        if debug:
            print ('## Debug: rep_list: (let_out, type_out)')
            print ('rep_list i:', i)
            print ('let_out', let_out)
            print ('type_out', type_out)
            
        out_list.append(let_out)
        out_list.append(type_out)

    return out_list
