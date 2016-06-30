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
import logging
import logging.handlers
import urlparse
import os
import sys
from datetime import datetime, timedelta
from pprint import pformat
import numpy as np
from pyresample import utils
from pyorbital import astronomy
from trollsched.spherical import get_twilight_poly
from trollsched.graph import Graph
from trollsched.satpass import get_next_passes, SimplePass
from trollsched.boundary import AreaDefBoundary
from trollsched.schedule import get_passes_from_xml_file

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
        pl.append(passes[s])
    
    n_vertices=1
    for g in grl:
        n_vertices *= g.order
    newgraph = Graph(n_vertices=n_vertices)

    print "newgraph order:", newgraph.order 
    
    newpasses = [tuple((0,None) for _ in range(len(statlst)))]
    parlist = [newpasses[0]]
    x=0
    while len(parlist) and x<5:
#         x+=1
        
#         print "-----------------------------------------"
#         print "with parlist",parlist
        
        newparlist = []
        for parnode in parlist:
            
            print "\nfor parnode",parnode,"from",parlist
            if parnode == ((None,None),(None,None)):
                print "skip parnode",parnode
                continue
            
            collected_newnodes = collect_nodes(0, parnode, grl, newgraph, newpasses, pl, delay)
            
            print "collected nodes",collected_newnodes
            
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

                wl = []
                s = 0
                for p,n in zip(parnode, newnode):
                    try:
                        
                        print "p,n",p,n
                        
                        if n[0] is None:
                            wl.append(0)
                        else:
                            if n[1] is None:
                                wl.append(grl[s].weight(p[0], n[0]))
                            else:
                                wl.append(n[1])
                        s += 1
                    except:
                        print "\nCATCH\nparnode",parnode,p,"newnode",newnode,n
                        print "\n",sys.exc_info()
                        raise
                
                we = sum(wl) # g.weight(p[0], n[0]) for n in newnode for p in parnode)
                ws = max(Counter(n[0] for n in parnode).values()) + max(Counter(n[0] for n in newnode).values()) - 2
                w = we / 2 ** ws
                
#                 print "weight:",parnode,newnode,"we",we,"ws",ws,"w",w

                if parnode == newpasses[0]:
                    
                    # for the starting point
#                     print "add_arc",0, newpasses.index(newnode), w
                    
                    newgraph.add_arc(0, newpasses.index(newnode) + 1, w)
                else:
                    
#                     print "add_arc",newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w
                    
                    newgraph.add_arc(newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w)

        parlist = newparlist
    
    print "leave add_graphs()"
    
    return statlst, newgraph, newpasses
    

# oder als Generator (mit yield statt return und for über parents statt rekursion)
def collect_nodes(statnr, parnode, graph_set, newgraph, newpasses, allpasses, delay=timedelta(seconds=0)):
    bufflist=[]
    p = parnode[statnr]
    g = graph_set[statnr]
    
    if p == (None,None):
        # there won't be any collectable nodes.
#         return [p]
        gn = [None]

    else:
    
        gn = g.neighbours(p[0])
    
    #     print "station",statnr,"parnode",p,"neighbours",gn
        
        if gn[0] > len(allpasses[statnr]):
            gn = [None]
            
    #         print "ENDE",statnr,p
        
    if  statnr + 1 == len(parnode):
        
        if p[1] is not None:
            
            bufflist.append( (p[0], None) )
            
#             print "simulated single ",p

        else:
#             print "single parnode",g.neighbours(p[0])
            bufflist = zip(gn, (None for _ in gn))
        
    else:
        
#         print "recurse"
        
        col = collect_nodes(statnr + 1, parnode, graph_set, newgraph, newpasses, allpasses)
        
#         print "returned:",col,"going into",gn

        if p[1] is not None:
            
#             print "simulating",p
            
            for x in col:
                bufflist.append( ((p[0], None), x) )
                
#                 print "simulated",(p, x)

        else:
            for n in gn:
            
#                 print "append",n,"+",col,"=>",
                for x in col: # zip(range(1,len(col)+1), col):
                    try:
                        
#                         if n==None:
#                             print "compare [",statnr,"][n",n,"] <=> [",1,"][x",x,"] none none",
#                         elif x == (None,None):
#                             print "compare [",statnr,"][n",n,"] <=> [",1,"][x",x,"]",allpasses[statnr][n-1],"none",
#                         else:
#                             print "compare [",statnr,"][n",n,"] <=> [",1,"][x",x,"]",allpasses[statnr][n-1],allpasses[statnr+1][x[0]-1],
                        
                        if n is None or x == (None,None) or allpasses[statnr][n-1].overlaps(allpasses[statnr+1][x[0]-1], delay):
#                             print "--> n=x"
                            bufflist.append( ((n, None), x) )
                            
                        elif allpasses[statnr][n-1].risetime > allpasses[statnr+1][x[0]-1].falltime + delay:
                            
#                             print "--> n>x","p",p,"n",n,"x",x,"w",g.weight(p[0], n)
                            
                            bufflist.append( ((n, g.weight(p[0], n)), x) )
                             
                        elif allpasses[statnr][n-1].falltime + delay < allpasses[statnr+1][x[0]-1].risetime:
                            
#                             print "--> n<x","p",parnode[statnr+1],"n",n,"x", x, "w",graph_set[statnr+1].weight(parnode[statnr+1][0], x[0])
                            
                            bufflist.append( ((n, None), (x[0], graph_set[statnr+1].weight(parnode[statnr+1][0], x[0]))) )
                            
#                         print bufflist[-1],
                    except:
                        print "\nCATCH\ngn:",gn,"-> n",n," col:",col,"-> x",x,"statnr",statnr,"statnr+i",statnr+1
                        print "len(allpasses -n -x)",len(allpasses[statnr]),len(allpasses[statnr+1])
                        print "allpasses-n",allpasses[statnr]
                        print "allpasses-x",allpasses[statnr+1]
                        print sys.exc_info()
                        raise
#                 print
#     print "returning",bufflist
    
    return bufflist



def get_combined_sched(allgraphs, allpasses):
    
    statlst, newgraph, newpasses = add_graphs(allgraphs, allpasses, delay=timedelta(seconds=300))

#     print_matrix(newgraph.adj_matrix, ly=5)
#     print_matrix(newgraph.weight_matrix, ly=5)
#     print newpasses

    dist, path = newgraph.dag_longest_path(0, len(newpasses))
    
    print "dist",dist,"path",path

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
    try:
        from trollsched.schedule import read_config
        import argparse
        #global logger
        logger = logging.getLogger("trollsched")
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--config", help="configuration file to use", default=None)
        parser.add_argument("-g", "--graph", help="save graph info to this directory", default=None)
        parser.add_argument("-r", "--report", help="xml report file", default=None)
        opts = parser.parse_args()

        if opts.config is None:
            parser.error("Configuration file required.")
        if opts.graph is None:
            parser.error("Graph directory required.")
        if opts.report is None:
            parser.error("Report file required.")
        
        # [coords, station, area, scores], forward, start
        station_list, forward, start = read_config(opts.config)
        graph = {}
        allpasses = {}
        for coords, station, area, scores in station_list:
            print "---",station,"---"
            graph[station] = Graph()
            graph[station].load(os.path.join(opts.graph, "graph." + station + ".npz"))
            
            print_matrix(graph[station].adj_matrix, ly=5)
            print_matrix(graph[station].weight_matrix, ly=5, lx=-1)
            
            allpasses[station] = get_passes_from_xml_file(os.path.join(opts.report, "acquisition-schedule-report." + station + ".xml"))
            print len(allpasses[station]),allpasses[station]
            
#             for v in graph[station].neighbours(1):
#                 print v, " : ", allpasses[station][v].risetime, "->", graph[station].weight(1, v)
        
        for s,g in graph.items():
            print "test folding",s
            test_folding(g)
        

        
        stats,schedule, (newgraph, labels) = get_combined_sched(graph, allpasses)

        print "stats",stats,"schedule",schedule,"labels",labels

        logger.debug(pformat(schedule))
        
        for opass in schedule:
            
            print opass,"->"
            
            for i,ipass in zip(range(len(opass)),opass):
                if ipass[0] is None:
                    continue
                print "->",i,ipass,allpasses[stats[i]][ipass[0]-1]
                allpasses[stats[i]][ipass[0]-1].rec = True
        
        logger.info("generating file")
    
    
        """
        TODO: 
        dateien für jede antenne erstellen.
        - scisys-schedule
        - xml-schedule
        funktionen ggf umschreiben/alternative impl
        """
    
#         if opts.scisys:
#             generate_sch_file(opts.scisys, station, allpasses[station], coords)

#         xmlfile = generate_xml_file(newpasses, start_time + timedelta(hours=start),
#                                     start_time + timedelta(hours=forward),
#                                     directory, station,
#                                     opts.report)


        
#         print_matrix(newgraph.adj_matrix, ly=5)
#         print_matrix(newgraph.weight_matrix, ly=5)

        print "test folding newgraph"
        test_folding(newgraph)
        
        print allpasses
        
        newgraph.save(os.path.join(opts.graph, "newgraph.comb"))
#         newgraph.export(labels=[str(label) for label in labels],
#                      filename=os.path.join(opts.graph,
#                                            "newsched.comb.gv"))

    except:
        logger.exception("Something wrong happened!")
        raise


