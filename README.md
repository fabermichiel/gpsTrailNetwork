gpsTrailNetwork
===============

creates a network of trails from several overlapping gps-tracks

This collection of scripts tries to create a network of trails (real world roads) from input gps-tracks.
Whenever one records different walks with a gps-receiver on the same path, the trails won't be the same. 
These scripts tries to average those trails and it tries to identify intersections. 
Output wil be a gpx-file with intersections as way- or routepoints and the connecting trails as routesegments
