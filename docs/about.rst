.. about_

About PyTroll-Schedule
======================

Case of One Receiving Station
-----------------------------
In the case of a single station, the procedure of scheduling is quite
straightforward. However, let us describe it in detail here, such that the
background will be set for the more complex case of multiple reception station
reception scheduling.

The first step to compute the schedule, is to know which satellites of interest
are going to be rising above the horizon during the duration of the schedule.
In order to find such cases, we retrieve the orbital information for each
satellite of interest and apply orbit prediction using the aiaa sgp4 algorithm
(ref). In practice, we use norad tle files (ref) as orbital elements, and the
python implementation of the sgp4 algorithm provided in pyorbital (ref). From
this, we then obtain a list of the coming overpasses for the station. We define
an overpass as a risetime and fall time for a given satellite, during which it
will be within reception reach of the station.

Now, we have to find the possible schedules for the station. The set of all
overpasses gives us all the reception possibilities for the station. However,
many of them will be in conflict with at least one other overpass and will be
a concurrent to the reception race. We say that two overpasses conflict when
the risetime dog one of them is comprised within the view time of the second.
In case of conflicts, the scheduling algorithm has to choose one or the other
overpass. However, in the case of several overpasses conflicting sequentially,
we have to find the possible paths through the conflicting zone. In order to do
that, we will use graph theory algorithms.

We define the graph of the conflicting zone with overpasses as vertices and
create an edge between two conflicting overpasses. To find the possible
non-conflicting combinations in this graph is actually searching for maximal
cliques in the complementary graph, for which we use the Bron-Kerbosch
algorithm.
#illustration click

we obtain thus groups of passes that are not conflicting in the time frame.
The next step is to find the optimal list of non conflicting passes under the
duration on the schedule.

Cases of Connected Stations
---------------------------
There are several ways to compute schedules for connected stations, two are
implemented in this program.

Several points should be considered:
* Technical equipment, reception of L-band, Ku-band, X-band?
* Geographic location, nearby or large distance between?

"Master-Slave" Operation
************************
The mode of master-slave is best suited for two stations, located next to each
other, with similar technical systems.

In this case a schedule for one, namely the "master" station, would be computed,
as if it were only this one station.

In a second step this schedule plan is used as a subtraction list when
computing the schedule for the second, the "slave" station.

Co-operating Stations
*********************
A mode of co-operating stations can consider the distance between different
geographical locations and differences in technical equipment, most notable
different reception capabilities (X- & L-band vs. L-band).

In this case, each station defines a time span requirement for each pass. Then,
if a connected station can fulfil this requirement and is scheduling the same
pass, we can say that the stations are redundant. To avoid such redundancy, we
can define ways to synchronise the schedule to optimise the intake of data and
fulfil the pareto condition.
A simple protocol can be used to perform this: both A and B provide alternatives
and compute the enhanced score for the schedule including the others pass.

B can delegate the pass only if it can assure that the time span requirement of
A is respected.

This operation can be extended to more than two stations, all receiving a
single-operation schedule and an individual cooperating-schedule.
