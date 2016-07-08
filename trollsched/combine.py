#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

# Author(s):

#   Alexander Maul <alexander.maul@dwd.de>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Combine several graphs.
"""
from datetime import datetime, timedelta
from trollsched.graph import Graph
from collections import Counter

def add_graphs(graphs, passes, delay=timedelta(seconds=0)):
    """add all graphs to one combined graph.
    """
    statlst = graphs.keys()

    print "station list"
    for s,g in graphs.items():
        print s,g.order
    
    
    grl = []
    pl = []
    for s in statlst:
        grl.append(graphs[s])
        pl.append(sorted(passes[s], key=lambda x: x.risetime)) # passes[s]) # 
    
    n_vertices=1
    for g in grl:
        n_vertices *= g.order
    newgraph = Graph(n_vertices=n_vertices)

    print "newgraph order:", newgraph.order 
    
    stopper = tuple( (None,None) for s in range(len(statlst)) )
    newpasses = [tuple( (pl[s][p-1],None) for s in range(len(statlst)) for p in grl[s].neighbours(0) )]
    parlist = [newpasses[0]]
    
    x=0
    while len(parlist) and x<5:
#         x+=1
        
#         print "-----------------------------------------"
#         print "with parlist",parlist
        
        newparlist = []
        for parnode in parlist:
            
#             print "\n---for parnode",parnode,"\n---from",parlist

            if parnode == stopper:
                # All antennas reached end of passes list in this line of 
                # possibilities.
                # stopper == ((None,None) * stations)
                #
                # If this happens for all elements of parlist, newparlist will 
                # stay empty and (at the bottom of this loop) replace parlist,
                # which as an empty list will cause the surrounding while-loop 
                # to end. 
                
#                 print "skip parnode eq stopper",parnode,stopper
                
                continue
            
            collected_newnodes = collect_nodes(0, parnode, grl, newgraph, newpasses, pl, delay)
            
#             print "collected nodes",collected_newnodes
            
            for newnode_list in collected_newnodes:
                newnode = tuple(newnode_list)
                if newnode not in newpasses:

#                     print "add newnode",newnode,"to newpasses",newpasses

                    newpasses.append(newnode)

#                     print "with index",newpasses.index(newnode)
#                     print "added",newnode,"index",newpasses.index(newnode)
                    
                if newnode not in newparlist:
                    newparlist.append(newnode)

#                     print "extended newparlist",newparlist
                
#                 print "collect weight:"
                
                # Collecting the weights from each stations weight-matrix ...
                wl = []
                s = 0
                for p,n in zip(parnode, newnode):
                    try:
                        
#                         print "p,n",p,n
                            
                        if n[0] is None:
                            wl.append(0)
                        else:
                            wl.append( n[1] or grl[s].weight(pl[s].index(p[0]) + 1, pl[s].index(n[0]) + 1) )
                        s += 1
                    except:
                        print "\nCATCH\nparnode",parnode,p,"\nnewnode",newnode,n
                        raise
                # sum of collected weights
                we = sum(wl)
                # vertices with reference to same sat-pass, could be 0, 1, 2.
                ws = 4 - count_neq_passes(parnode) - count_neq_passes(newnode)
                # apply vertix-count to weight-sum
                w = we / 2 ** ws
                
#                 print "weight:\n",parnode,"\n",newnode,"\nwl",wl,"we",we,"ws",ws,"w",w

                if parnode == newpasses[0]:
                    
                    # for the starting point
#                     print "add_arc",0, newpasses.index(parnode), w
                    
                    newgraph.add_arc(0, newpasses.index(parnode) + 1, w)
#                 else:
                    
#                     print "add_arc",newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w
                    
                newgraph.add_arc(newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w)

        parlist = newparlist
    
    print "leave add_graphs()"
    
    return statlst, newgraph, newpasses
    

def count_neq_passes(pl):
    if len(pl):
        r = []
        s = 1
        for q in pl[1:]:
            if pl[0] != q:
                r.append(q)
        s += count_neq_passes(r)
        return s
    else:
        return 0


# oder als Generator (mit yield statt return und for über parents statt rekursion)
def collect_nodes(statnr, parnode, graph_set, newgraph, newpasses, passes_list, delay=timedelta(seconds=0)):
    bufflist=[]
    p = parnode[statnr]
    g = graph_set[statnr]
    
    if p is (None,None):
        
#         print "there won't be any collectable nodes."
        
        gn = [None]

    else:
    
        gn = g.neighbours(passes_list[statnr].index(p[0]) + 1)
    
#         print "station",statnr,"parnode",p,"neighbours",gn
        
        if gn[0] > len(passes_list[statnr]):
            gn = [None]
            
#             print "ENDE",statnr,p
        
    if  statnr + 1 == len(parnode):
        # It's the 'rightmost' of the list parnode, 
        # and the deepest point of the recursion.
        if p[1] is not None:
            # If the parent node is a simulated vertix it is appended.
            bufflist.append( (p[0], None) )
            
#             print "simulated single ",p

        elif None in gn:
            # That's "no further connection".
            bufflist = [(None, None)]
            
        else:
            # prepare to return just the list of neighbouring vertices.

#             print "single parnode",gn
            
            bufflist = zip((passes_list[statnr][n-1] for n in gn), (None for _ in gn))
        
    else:
        # If it's not the last element of the list parnode, we recurse and then
        # permutade all vertix-lists together.

#         print "recurse"
        
        col = collect_nodes(statnr + 1, parnode, graph_set, newgraph, newpasses, passes_list)
        
#         print "returned:",col,"going into",gn

        if p[1] is not None:
            # A simulated parent node is permutated with all neighbours.
            
#             print "simulating",p
            
            for x in col:
                bufflist.append( ((p[0], None), x) )
                
#                 print "simulated",(p, x)

        else:
            # Creating the permutation of all neighbours with the list returned
            # by the recursion. 
            for n in gn:
            
#                 print "append",n,"+",col,"=>",
                for cx in col:
                    try:

#                         if n is None:
#                             print "compare [",statnr,"][n",n,"] <=>",
#                         else:
#                             print "compare [",statnr,"][n",passes_list[statnr][n-1],"] <=>",
#                         print "[",1,"][cx",cx,"]"
                        
                        if n is None:
                            # The end-of-neighbours dummy.
                            
#                             print "XXXXXXXXXXXXXXXXXXXXXXX"
                            
                            bufflist.append( ((None, None), cx) )
                            
                        elif cx == (None,None) or passes_list[statnr][n-1].overlaps(cx[0], delay):
                            # Two passes overlapping, no special handling required.

#                             print "--> n=cx"
                            
                            bufflist.append( ((passes_list[statnr][n-1], None), cx) )
                            
                        elif passes_list[statnr][n-1].risetime > cx[0].falltime + delay:
                            # If the current parent node's pass is not overlapping 
                            # but AFTER the pass from the recursion-return-list
                            # the current parent node gets "simulated".
                            
#                             print "--> n>cx","p",p,"n",n,"cx",cx,"w",g.weight(passes_list[statnr].index(p[0]) + 1, n)
                            
                            bufflist.append( ((passes_list[statnr][n-1], g.weight(passes_list[statnr].index(p[0]) + 1, n)), cx) )
                             
                        elif passes_list[statnr][n-1].falltime + delay < cx[0].risetime:
                            # If the current parent node's pass is not overlapping 
                            # but BEFORE the pass from the recursion-return-list
                            # the recursion-list-node gets "simulated".
                            
#                             print "--> n<cx","p",parnode[statnr+1],"n",n,"cx", cx, "w",graph_set[statnr+1].weight(
#                                                         passes_list[statnr+1].index(parnode[statnr+1][0]), 
#                                                         passes_list[statnr+1].index(cx[0])
#                                                         )
                            
                            bufflist.append( (
                                        (passes_list[statnr][n-1], None), 
                                        (cx[0], graph_set[statnr+1].weight(
                                                    passes_list[statnr+1].index(parnode[statnr+1][0]) + 1, 
                                                    passes_list[statnr+1].index(cx[0]) + 1
                                                    )) 
                                        ) )
                        
                        else:
                            print "XXXXXXXXXXXXXXXXXXXXXXX something curious happened ..."
                        
#                         print bufflist[-1],
                    except:
                        print "\nCATCH\ngn:",gn,"-> n",n," col:",col,"-> cx",cx,"statnr",statnr,"statnr+i",statnr+1
                        print "len(passes_list -n -cx)",len(passes_list[statnr]),len(passes_list[statnr+1])
                        print "passes_list-n",passes_list[statnr]
                        print "passes_list-x",passes_list[statnr+1]
                        raise
#                 print
#     print "returning",bufflist
    
    return bufflist



def get_combined_sched(allgraphs, allpasses, delay_sec=60):
    
    delay = timedelta(seconds=delay_sec)

#     for s,g in allgraphs.items():
#         print "***",s,"***"
#         print_matrix(g.adj_matrix, ly=5)
#         print_matrix(g.weight_matrix, ly=5, lx=-1)
    
    statlst, newgraph, newpasses = add_graphs(allgraphs, allpasses, delay)

#     print_matrix(newgraph.adj_matrix, ly=5)
#     print_matrix(newgraph.weight_matrix, ly=5)
#     print newpasses

    dist, path = newgraph.dag_longest_path(0, len(newpasses))
    
    print "---dist---",dist
    print "---path---",path

    del dist
    return statlst,[newpasses[idx - 1] for idx in path[1:-1]], (newgraph, newpasses)






def print_matrix(m, ly = -1, lx = -1):
    for i,l in zip(range(ly), m[0:ly]):
        print i,":",l[:lx]
    print "[..., ...]"
    for i,l in zip(range(len(m)-ly-1, len(m)-1), m[-ly:]):
        print i,":",l[:lx]



def test_folding(g):
    for u in range(g.order):
        for n in g.neighbours(u):
            if n<u:
                print n,"<",u



if __name__ == '__main__':
    
    import logging
    import logging.handlers
    import urlparse
    import os
    from pprint import pformat
    import numpy as np
    from pyresample import utils
    from pyorbital import astronomy
    from trollsched.spherical import get_twilight_poly
    from trollsched.satpass import get_next_passes, SimplePass
    from trollsched.boundary import AreaDefBoundary
    from trollsched.schedule import get_passes_from_xml_file,generate_sch_file,generate_xml_file,parse_datetime
    from trollsched.schedule import combined_stations

    try:
        from trollsched.schedule import read_config
        import argparse
        #global logger
        logger = logging.getLogger("trollsched")
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--config", help="configuration file to use", default=None)
        parser.add_argument("-g", "--graph", help="save graph info to this directory", default=None)
        parser.add_argument("-r", "--report", help="xml report file", default=None)
        parser.add_argument("-s", "--start-time",
                            type=parse_datetime,
                            help="start time of the schedule to compute")
        parser.add_argument("--scisys", help="scisys schedule file", default=None)
        opts = parser.parse_args()

        if opts.config is None:
            parser.error("Configuration file required.")
        if opts.graph is None:
            parser.error("Graph directory required.")
        if opts.report is None:
            parser.error("Report file required.")
        if opts.start_time:
            start_time = opts.start_time
        else:
            start_time = datetime.utcnow()
        
        # [coords, station, area, scores], forward, start
        station_list, forward, start = read_config(opts.config)
        graph = {}
        allpasses = {}
        for coords, station, area, scores in station_list:
            graph[station] = Graph()
            graph[station].load(os.path.join(opts.graph, "graph." + station + ".npz"))

#             print "---",station,"---"
#             print_matrix(graph[station].adj_matrix, ly=5)
#             print_matrix(graph[station].weight_matrix, ly=5, lx=-1)
            
#             allpasses[station] = get_passes_from_xml_file(os.path.join(opts.report, "acquisition-schedule-report." + station + ".xml"))
#             print len(allpasses[station]),allpasses[station]
            
#             for v in graph[station].neighbours(1):
#                 print v, " : ", allpasses[station][v].risetime, "->", graph[station].weight(1, v)
        
        import pickle
        ph = open(os.path.join(opts.graph,"opts.pkl"), "rb")
        opts = pickle.load(ph)
        ph.close()
        ph = open(os.path.join(opts.graph,"allpasses.pkl"), "rb")
        allpasses = pickle.load(ph)
        ph.close()
        
        
        combined_stations(opts, station_list,graph, allpasses, start_time, start, forward)

        
    except:
        logger.exception("Something wrong happened!")
        raise


