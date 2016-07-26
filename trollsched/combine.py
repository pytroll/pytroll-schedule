#!/usr/bin/env python
# -*- coding: utf-8 -*-
from scipy.constants.constants import year

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

def add_graphs(graphs, passes, delay=timedelta(seconds=0)):
    """add all graphs to one combined graph.
    """
    statlst = graphs.keys()

    def count_neq_passes(pl):
        """Counts how many satellite passes in a list are really distinct (satellite/epoch).

        TODO: The "same epoch" is only guessed by comparing with a time window
        hard-coded in SimplePass.__eq__() -- this is highly impovable!
        """
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

    print "station list"
    for s, g in graphs.items():
        print s, g.order

    # graphs and allpasses are hashmaps of sets or so, but we need lists of lists,
    # forthat they are copied.
    grl = []
    pl = []
    for s in statlst:
        grl.append(graphs[s])
        pl.append(sorted(passes[s], key=lambda x: x.risetime))

    # rough estimate for the size of the combined passes' graph.
    n_vertices = 1
    for g in grl:
        n_vertices += g.order
    n_vertices *= len(statlst)
    newgraph = Graph(n_vertices=n_vertices)

    print "newgraph order:", newgraph.order

    # this value signals the end, when no more passes from any antenna are available.
    stopper = tuple((None, None) for s in range(len(statlst)))

    # the new passes list, it'll be filled with tuples of one pass per antenna.
    # it's initialized with the first passes.
    newpasses = [tuple((pl[s][p - 1], None) for s in range(len(statlst)) for p in grl[s].neighbours(0))]
    parlist = [newpasses[0]]

    x = 0
    while len(parlist) and x < 5:
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
                for s, p, n in zip(range(len(statlst)), parnode, newnode):
                    try:
                        if n[0] is None:
                            wl.append(0)
                        else:
                            wl.append(n[1] or grl[s].weight(pl[s].index(p[0]) + 1, pl[s].index(n[0]) + 1))
                    except:
                        print "\nCATCH\nstat", s, "parnode", parnode, p, "\nnewnode", newnode, n
                        raise
                # sum of collected weights
                ws = sum(wl)
                # vertices with reference to same sat-pass, 'we' can result to 0, 1, 2.
                we = 4 - count_neq_passes(parnode) - count_neq_passes(newnode)
                # apply vertix-count to weight-sum
                w = ws / 2 ** we

#                 print "weight:\n",parnode,"\n",newnode,"\nwl",wl,"-> ws",ws,"/2**we",we,"== w",w

                if parnode == newpasses[0]:

                    # for the starting point
#                     print "add_arc",0, newpasses.index(parnode), w

                    newgraph.add_arc(0, newpasses.index(parnode) + 1, w)
#                 else:

#                     print "add_arc",newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w

                newgraph.add_arc(newpasses.index(parnode) + 1, newpasses.index(newnode) + 1, w)

        parlist = newparlist

    print "len(newpasses) =", len(newpasses)
    print "leave add_graphs()"

    return statlst, newgraph, newpasses


def collect_nodes(statnr, parnode, graph_set, newgraph, newpasses, passes_list, delay=timedelta(seconds=0)):
    """Collect all nodes reachable from the nodes in parnode, creating all combinations.

    RETURN: [[a1, b1], [a1, b2], ..., [a2, b1], ...]
    """
    bufflist = []
    p = parnode[statnr]
    g = graph_set[statnr]

    def overlap_any(this, test_list):
        """Tests if this overlapps any of the new-nodes in test_list.

        The new-nodes are in form (vertix, simulated-weight), only nodes without
        simulated weight are considered in the test.

        RETURN: -1 | 0 | +1 , if this lies before, overlapps any, or lies after
        the nodes in test_list.
        """
        if len(test_list) == 0:
            return 0
        minrise = datetime.utcnow() + timedelta(days=300)
        maxfall = datetime.utcnow() - timedelta(days=300)
        for p in test_list:
            if p[0] is None or p[1] is not None:
                continue
            if p[0].risetime < minrise:
                minrise = p[0].risetime
            if p[0].falltime > maxfall:
                maxfall = p[0].falltime
        if minrise > maxfall:
            return 0
        elif this.falltime < minrise:
            return -1
        elif this.risetime > maxfall:
            return +1
        else:
            return 0

    if p == (None, None):
        # There won't be any collectable nodes.
        # This None will act as a filler in the combined-vertices-tuples,
        # to get the access-by-index right.

#         print "there won't be any collectable nodes."

        gn = [None]

    elif p[1] is not None:
        # A simulated parent node is set as neighbours' list.
        # It'll be processed as if it's the node which occurs in this
        # time-slot -- which it propably does, otherwise it's subjected
        # to simulation (again!).

        gn = [passes_list[statnr].index(p[0]) + 1]

#         print "simulated",p, gn

    else:

        gn = g.neighbours(passes_list[statnr].index(p[0]) + 1)

#         print "station",statnr,"parnode",p,"neighbours",gn

        if gn[0] > len(passes_list[statnr]):
            gn = [None]

#             print "ENDE",statnr,p

#     print "statnr",statnr,"len(parnode)",len(parnode),parnode,"[stat]->",p

    if  statnr + 1 == len(parnode):
        # It's the 'rightmost' of the list parnode,
        # and the deepest point of the recursion.

        if None in gn:
            # That's "no further connection".
            # It get's a special treatment, because there is no None in the
            # passes-list we could access by index.

            bufflist = [[(None, None)]]

        else:
            # Prepare to return just the list of neighbouring vertices.

#             print "single parnode",gn

            for m in zip((passes_list[statnr][n - 1] for n in gn), (None for _ in gn)):
                bufflist.append([m])

    else:
        # Since it's not the last element of the list parnode, we recurse and
        # then permutade all vertix-lists together.

#         print "recurse"

        col = collect_nodes(statnr + 1, parnode, graph_set, newgraph, newpasses, passes_list)

#         print "returned:",col,"going into",gn

        # Creating the permutation of all neighbours with the list returned
        # by the recursion.
        # A simulated parent node is seen as a regular list of neighbours.
        for n in gn:

#                 print "append",n,"+",col,"=>"

            for cx in col:
                try:

                    if n is None:
                        # The end-of-neighbours dummy.

#                             print "XXXXXXXXXXXXXXXXXXXXXXX"

                        cc = cx[:]
                        cc.insert(0, (None, None))
                        bufflist.append(cc)

                    else:

#                         if n is None:
#                             print "compare [",statnr,"][n",n,"] <=>",
#                         else:
#                             print "compare [",statnr,"][n",passes_list[statnr][n-1],"] <=>",
#                         print "[",1,"][cx",cx,"]"

                        overlap = overlap_any(passes_list[statnr][n - 1], cx)

                        if overlap == 0: # cx[0] == (None, None) or
                            # Two passes overlapping, no special handling required.

#                             print "--> n=cx"

                            cc = cx[:]
                            cc.insert(0, (
                                          (passes_list[statnr][n - 1], None)
                                          ))
                            bufflist.append(cc)

                        elif overlap > 0:
                            # If the current parent node's pass is not overlapping
                            # but AFTER the pass from the recursion-return-list
                            # the current parent node gets "simulated".

#                             print "--> n>cx","p",p,"n",n,"cx",cx,"w",g.weight(passes_list[statnr].index(p[0]) + 1, n)

                            cc = cx[:]
                            cc.insert(0, (
                                          passes_list[statnr][n - 1],
                                          g.weight(passes_list[statnr].index(p[0]) + 1, n)
                                          ))
                            bufflist.append(cc)

                        elif overlap < 0:
                            # If the current parent node's pass is not overlapping
                            # but BEFORE the pass from the recursion-return-list
                            # the recursion-list-node gets "simulated".
                            cc = [
                                  (c[0], graph_set[s].weight(
                                            passes_list[s].index(parnode[s][0]) + 1,
                                            passes_list[s].index(c[0]) + 1
                                            )
                                  ) if c != (None, None) else (None, None)
                                  for s, c in zip(range(statnr + 1, len(parnode)), cx)
                                ]
                            cc.insert(0, (passes_list[statnr][n - 1], None))

#                             print "--> n<cx","p",passes_list[statnr][n-1],"n",n,"cx", cx, "cc",cc

                            bufflist.append(cc)

                        else:
                            print "uh-oh, something curious happened ..."

    #                         print bufflist[-1],
                except:
                    print "\nCATCH\ngn:", gn, "-> n", n, " col:", col, "-> cx", cx, "statnr", statnr, "statnr+i", statnr + 1
                    print "len(passes_list -n -cx)", len(passes_list[statnr]), len(passes_list[statnr + 1])
                    for s in range(statnr, len(passes_list)):
                        print "passes_list[", s, "] =>", passes_list[s]
                    raise

#                 print
#     print "returning",bufflist
#     print

    return bufflist



def get_combined_sched(allgraphs, allpasses, delay_sec=60):

    delay = timedelta(seconds=delay_sec)

#     for s, g in allgraphs.items():
#         print "***", s, "***"
#         print_matrix(g.adj_matrix, ly=5)
#         print_matrix(g.weight_matrix, ly=5, lx=-1)

    statlst, newgraph, newpasses = add_graphs(allgraphs, allpasses, delay)

#     print_matrix(newgraph.adj_matrix, ly=5)
#     print_matrix(newgraph.weight_matrix, ly=5)
#     print newpasses

#     for s, g in allgraphs.items():
#         print "test folding", s
#         test_folding(g)
#     print "test folding newgraph"
#     if test_folding(newgraph):
#         print_matrix(newgraph.adj_matrix, 25, 36)

    dist, path = newgraph.dag_longest_path(0, len(newpasses))

    print "---dist---", dist
    print "---path---", path

    del dist
    return statlst, [newpasses[idx - 1] for idx in path[1:-1]], (newgraph, newpasses)






def print_matrix(m, ly=-1, lx=-1):
    """For DEBUG: Prints one of the graphs' backing matrix without
    flooding the screen.

    It'll print the first lx columns from the first ly rows, then
    the last lx columns from the last ly rows.
    """
    for i, l in zip(range(ly), m[0:ly]):
        print i, ":", l[:lx], "..."
    print "[..., ...]"
    for i, l in zip(range(len(m) - ly - 1, len(m) - 1), m[-ly:]):
        print i, ": ...", l[-lx:]



def test_folding(g):
    r = False
    for u in range(g.order):
        for n in g.neighbours(u):
            if n < u:
                print n, "<", u
                r = True
    return r



if __name__ == '__main__':

    import logging
    import logging.handlers
    import urlparse
    import os
    import pickle
    from pprint import pformat
    import numpy as np
    from pyresample import utils
    from pyorbital import astronomy
    from trollsched.spherical import get_twilight_poly
    from trollsched.satpass import get_next_passes, SimplePass
    from trollsched.boundary import AreaDefBoundary
    from trollsched.schedule import get_passes_from_xml_file, generate_sch_file, generate_xml_file, parse_datetime
    from trollsched.schedule import combined_stations, build_filename

    try:
        from trollsched.schedule import read_config
        import argparse
        # global logger
        logger = logging.getLogger("trollsched")
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--config", default=None,
                            help="configuration file to use")
        parser.add_argument("-s", "--start-time", type=parse_datetime,
                            help="start time of the schedule to compute")
        parser.add_argument("-o", "--output-dir", default=None,
                            help="where to put generated files")
        parser.add_argument("-x", "--xml", action="store_true",
                           help="generate an xml request file (schedule)"
                           )
        parser.add_argument("-r", "--report", action="store_true",
                           help="generate an xml report file (schedule)")
        parser.add_argument("--scisys", action="store_true",
                           help="generate a SCISYS schedule file")
        parser.add_argument("-p", "--plot", action="store_true",
                            help="generate plot images")
        parser.add_argument("-g", "--graph", action="store_true",
                            help="save graph info")
        opts = parser.parse_args()

        if opts.config is None:
            parser.error("Configuration file required.")
        if opts.start_time:
            start_time = opts.start_time
        else:
            start_time = datetime.utcnow()

        # [coords, station, area, scores], forward, start, pattern
        station_list, forward, start, pattern = read_config(opts.config)

        pattern_args = {
                        "output_dir":opts.output_dir,
                        "date":start_time.strftime("%Y%m%d"),
                        "time":start_time.strftime("%H%M%S")
                        }
        dir_output = build_filename("dir_output", pattern, pattern_args)
        if not os.path.exists(dir_output):
            print dir_output,"does not exist!"
            sys.exit(1)
        ph = open(os.path.join(dir_output, "opts.pkl"), "rb")
        opts = pickle.load(ph)
        ph.close()

        graph = {}
        allpasses = {}
        for coords, station, area, scores in station_list:
            graph[station] = Graph()
            graph[station].load(build_filename("file_graph", pattern, pattern_args))

#             print "---",station,"---"
#             print_matrix(graph[station].adj_matrix, ly=5)
#             print_matrix(graph[station].weight_matrix, ly=5, lx=-1)

#             allpasses[station] = get_passes_from_xml_file(os.path.join(opts.report, "acquisition-schedule-report." + station + ".xml"))
#             print len(allpasses[station]),allpasses[station]

#             for v in graph[station].neighbours(1):
#                 print v, " : ", allpasses[station][v].risetime, "->", graph[station].weight(1, v)

            ph = open(os.path.join(build_filename("dir_output", pattern, pattern_args), "allpasses.%s.pkl" % station), "rb")
            allpasses[station] = pickle.load(ph)
            ph.close()


        from trollsched.schedule import conflicting_passes
        totpas = []
        for s, sp in allpasses.items():
            print "len(sp)", s, len(sp)
            totpas.extend(list(sp))
        passes = sorted(totpas, key=lambda x: x.risetime)
        cpg = conflicting_passes(passes, timedelta(seconds=600))
        print "ALLPASSES", len(allpasses) # ,allpasses
        print "PASSES", len(passes) # ,passes
        print "CONFLGRPS", len(cpg) # ,cpg
        print "MAX", max([len(g) for g in cpg])



        combined_stations(opts, station_list, graph, allpasses, start_time, start, forward)


    except:
        logger.exception("Something wrong happened!")
        raise


