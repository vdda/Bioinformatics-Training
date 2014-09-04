#!/usr/bin/env python
"""
A simple script to convert Celera(R) Assembler's "best.edges" to a GEXF which can be used to feed into Gephi to
check the topology of the best overlapping graph.

help:
    python CA_best_edge_to_GML.py -h

Require networkx python module and Celera Assembler 8.1 in the path, defaults assume Celera Assembler was ran as part of
a SMRT Analysis workflow 
"""

import networkx as nx
import os
import shlex
import sys, argparse
import subprocess

def main(argv):

    parser = argparse.ArgumentParser(description='convert Celera(R) Assembler\'s \"best.edges\" to a gexf graph file')

    parser.add_argument('-g','--gkp_store', help='CA gkp_store directory, (celera-assembler.gkpStore)', default="celera-assembler.gkpStore")
    parser.add_argument('-t','--tig_store', help='CA tig_store directory, (celera-assembler.tigStore)', default="celera-assembler.tigStore")
    parser.add_argument('-b','--best_edge', help='CA best edge file, (./4-unitigger/best.edges)', default="./4-unitigger/best.edges")
    parser.add_argument('-c','--csv_data',  help='file containing arbitrary data in csv format', required=False)
    parser.add_argument('-o','--output',    help='output gexf file, (output.gexf)', default="output")

    args = parser.parse_args()

    gkp_store = args.gkp_store
    tig_store = args.tig_store
    best_edge = args.best_edge
    csv = args.csv_data
    output = args.output

    G=nx.DiGraph()
    frg_to_tig = {}
    cout = {}
    args = shlex.split("tigStore -g %s -t %s 2 -D unitiglist" % (gkp_store, tig_store ))
    out = subprocess.check_output(args)
    out = out.split("\n")
    for l in out:
        l = l.strip().split()
        if len(l) == 0: continue
        if l[0] == "maID": continue
        unitig_id = int(l[0])

        os.system("tigStore -g %s -t %s 2 -d frags -u %d > frag_list" % ( gkp_store, tig_store, unitig_id) )

        args = shlex.split( "tigStore -g %s -t %s 2 -d frags -u %d" % ( gkp_store, tig_store, unitig_id) )
        f_out = subprocess.check_output(args)
        f_out = f_out.split("\n")
        for l in f_out:
            """FRG    1453 179419,182165"""
            l = l.replace(",", " ")
            l = l.strip().split()
            if len(l) == 0: continue
            frg_id = l[1]
            frg_to_tig[frg_id] = unitig_id
    if(csv):
        with open(csv) as fin:
            for l in fin:
                l = l.strip().split(",")
                contig, cov, size, ref = l
                cout[contig] = size, cov, ref


    with open(best_edge) as f:
        for l in f:
            if l[0] == "#": continue
            l = l.strip().split()
            id1, lib_id, best5, o1, best3, o3, j1, j2 = l
    #        id1, lib_id, best5, o1, best3, o3 = l
            try:
                G.add_node(id1, label="utg%s" % frg_to_tig[id1], size=int(cout["unitig_%s"%frg_to_tig[id1]][0]), cov=float(cout["unitig_%s"%frg_to_tig[id1]][1]),ref=(cout["unitig_%s"%frg_to_tig[id1]][2]))
            except KeyError:
                G.add_node(id1, label="utg%s" % frg_to_tig[id1], size=int(0), cov=float(0))
            if best5 != "0":
                G.add_edge(best5, id1)
            if best3 != "0":
                G.add_edge(id1, best3)

    output_gexf = "%s.gexf" % output
    nx.write_gexf(G, output_gexf)

if __name__ == "__main__":
   main(sys.argv[1:])
