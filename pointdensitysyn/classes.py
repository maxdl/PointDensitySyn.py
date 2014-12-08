import sys
import exceptions
import os.path
import random
import itertools
import unicode_csv
import geometry
from fileIO import *
from stringconv import *


# Convenience functions
def dotProgress(x, linelength=80, char='.'):
    """ Simple progress indicator on sys.stdout
    """
    sys.stdout.write(char)
    if (x + 1) % linelength == 0:
        sys.stdout.write('\n')

#
# Classes
#
class Point(geometry.Point):
    def __init__(self, x=None, y=None, ptype="", profile=None):
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile != None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.skipped = False
        self.type = ptype
        self.cluster = None

    def determineStuff(self):
        """ Calculates various properties of self; enables control of when
            time-consuming calculations of properties are performed during 
            program execution
        """
        if self.isWithinHole:
            if self.type == "point":
                ProfileMessage(self.profile, "Discarding point at %s: Located "
                                             "within a profile hole" % self)
            self.skipped = True
            return
        if not hasattr(self, "distToPath"):
            self.distToPath = self.getDistToPath()
        if not self.isWithinProfile:
            if not self.isWithinShell:
                if self.type == "point":
                    ProfileMessage(self.profile, "Discarding point at %s: "
                                            "Located outside the shell" % self)
                self.skipped = True
                return
        self.lateralDistPSD_AZ, self.normLateralDistPSD_AZ,\
                self.closestPSD_AZ = self.lateralDistPSD_AZ()
        self.straightDistToPSD_AZ = self.straightDistSyn(self.profile.path,
                                                         self.profile.totalSynm)
        if not (hasattr(self, "isWithinPSD_AZ") and
                hasattr(self, "isAssociatedWithPSD_AZ") and
                hasattr(self, "associatedPSD_AZ")):
            self.isWithinPSD_AZ, self.isAssociatedWithPSD_AZ,\
                    self.associatedPSD_AZ  = self.psd_az_association()


    @property
    def isWithinHole(self):
        """  Determine whether self is inside a profile hole
        """
        is_within_hole = False
        for h in self.profile.path.holeli:
            if self.isWithinPolygon(h):
                is_within_hole = True
                break
        return is_within_hole

    @property
    def isWithinProfile(self):
        if (not self.isWithinPolygon(self.profile.path)) or self.isWithinHole:
            return False
        else:
            return True

    @property
    def getIsWithinPSD_AZ(self):
        if not hasattr(self, "isWithinPSD_AZ"):
            self.isWithinPSD_AZ, self.isAssociatedWithPSD_AZ,\
                    self.associatedPSD_AZ  = self.psd_az_association()
        return self.isWithinPSD_AZ

    @property
    def getIsAssociatedWithPSD_AZ(self):
        if not hasattr(self, "isAssociatedWithPSD_AZ"):
            self.isWithinPSD_AZ, self.isAssociatedWithPSD_AZ,\
                    self.associatedPSD_AZ  = self.psd_az_association()
        return self.isAssociatedWithPSD_AZ


    @property
    def isWithinShell(self):
        if not hasattr(self, "distToPath"):
            self.distToPath = self.getDistToPath()
        return (not self.isWithinProfile and
                 abs(self.distToPath) < self.profile.shell_width_in_pixels)

    @property
    def isAssociatedWithPath(self):
        if not hasattr(self, "distToPath"):
            self.distToPath = self.getDistToPath()
        return abs(self.distToPath) <= self.profile.spatial_resolution_in_pixels

    @property
    def isAssociatedWithProfile(self):
        return self.isWithinProfile or self.isAssociatedWithPath

    def getDistToPath(self):
        d = self.perpendDistClosedPath(self.profile.path)
        if not self.isWithinProfile:
            return -d if d > 0 else d
        else:
            return d if d > 0 else -d

    def lateralDistPSD_AZ(self):
        mindist = sys.maxint
        closestPSD_AZ = None
        for psdaz in self.profile.psd_az_li:
            d = (self.lateralDistSyn(self.profile.path, psdaz.synm) -
                    (psdaz.synm.length / 2))
            if d < mindist:
                mindist = d
                closestPSD_AZ = psdaz
        if not closestPSD_AZ:
            raise ProfileError(self.profile, "could not determine lateral "
                                "distance to a PSD/active zone of point at %s"
                                % self)
        mindist = self.lateralDistSyn(self.profile.path, closestPSD_AZ.synm)
        normdist = mindist / (closestPSD_AZ.synm.length / 2)
        return mindist, normdist, closestPSD_AZ

    def psd_az_association(self):
        isWithinPSD_AZ = False
        isAssociatedWithPSD_AZ = False
        associatedPSD_AZ = None
        mindist = sys.maxint
        for psd_az in self.profile.psd_az_li:
            if self.isWithinPolygon(psd_az.outline):
                isWithinPSD_AZ = True
                isAssociatedWithPSD_AZ = True
                associatedPSD_AZ = psd_az
                break
            dist = self.perpendDistClosedPath(psd_az.outline, 
                                              doNotCareIfOnOrOffSeg=True)
            if dist <= self.profile.spatial_resolution_in_pixels:
                isAssociatedWithPSD_AZ = True
                if dist < mindist:
                    associatedPSD_AZ = psd_az
                    mindist = dist
        return isWithinPSD_AZ, isAssociatedWithPSD_AZ, associatedPSD_AZ


    def determineNearestNeighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points are
        # in pointli
        self.nearestNeighbourDist = None
        mindist = float(sys.maxint)
        for p in pointli:
            # I will instead exclude non-desired points from the supplied point
            # list *before* calling this function
            #if p is not self and p.isAssociatedWithProfile:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxint):
            return None
        else:
            self.nearestNeighbourDist = mindist
            return self.nearestNeighbourDist
                
    def determineNearestLateralNeighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points are
        # in pointli
        self.nearestLateralNeighbourDist = None
        mindist = float(sys.maxint)
        for p in pointli:
            if p is not self:
                d = self.lateralDistToPoint(p, self.profile.path)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxint):
            return None
        else:
            self.nearestLateralNeighbourDist = mindist
            return self.nearestLateralNeighbourDist

        
    def lateralDistToPoint(self, p2, propath):
        """ Determine lateral distance to a point p2 along profile path.
        """
        path = geometry.SegmentedPath()
        p2_project, p2_seg_project = p2.projectOnClosedPath(propath)
        project, seg_project = self.projectOnClosedPath(propath)
        path.extend([project, p2_project])
        if p2_seg_project < seg_project:
            path.reverse()
        for n in range(min(p2_seg_project, seg_project) + 1,
                       max(p2_seg_project, seg_project)):
            path.insert(len(path)-1, propath[n])
        L = path.length
        return min(L, propath.perimeter - L)

    def lateralDistSyn(self, path, sm):
        """ Determine lateral distance to center of synaptic membrane sm, which
            is a subpath of path.
        """
        return self.lateralDistToPoint(sm.centerPoint, path)

    def straightDistSyn(self, pm, sm):
        """ Determine distance of self to synaptic membrane sm; if self's
            projection is not on sm, return distance to the nearest endnode
            of sm.
        """
        d = self.perpendDist(sm)
        if d == None:
            d = min(abs(self.dist(sm[0])), abs(self.dist(sm[-1])))
        if self.isWithinPolygon(pm):
            return d
        else:
            return -d


class PlasmaMembraneData(geometry.SegmentedPath):
    def __init__(self, pointlist=[], profile=None):
        geometry.SegmentedPath.__init__(self, pointlist)
        self.profile = profile
        self.holeli = []

        
    def addHole(self, pointlist=[]):
        self.holeli.append(pointlist)

    @property
    def area(self):
        """ Determine area of profile, excluding holes
        """
        return (super(PlasmaMembraneData, self).area -
                sum([h.area for h in self.holeli]))

        
    def contains(self, p):
        """  Determine whether point p is inside profile border 
        """
        if not p:
            return None
        return p.isWithinProfile


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        try:
            self.extend([Point(p.x, p.y, ptype, profile) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError, 'not a list of Point elements'

    def determineStuff(self):
        """ Determines properties of points in self; not really necessary, but
            used to avoid time-consuming calculations when saving output
        """
        for p in self:
            p.determineStuff()
        self[:] = [p for p in self if not p.skipped]


class ClusterData(list):
    def __init__(self, pli=[], profile=None):
        try:
            self.extend([Point(p.x, p.y, ptype="point", profile=profile)
                        for p in pli])
        except (AttributeError, IndexError):
            raise TypeError, 'not a Point list'
        self.profile = profile
        self.convexHull = geometry.SegmentedPath()

    def lateralDist(self, c2, pm):
        """ Determine lateral distance to a cluster c2 along plasma membrane.
        """
        path = geometry.SegmentedPath()
        c2_project, c2_seg_project = c2.convexHull.centroid.projectOnClosedPath(pm)
        project, seg_project = self.convexHull.centroid.projectOnClosedPath(pm)
        path.extend([project, c2_project])
        if c2_seg_project < seg_project:
            path.reverse()
        for n in range(min(c2_seg_project, seg_project) + 1,
                       max(c2_seg_project, seg_project)):
            path.insert(len(path)-1, pm[n])
        L = path.length
        return min(L, pm.perimeter - L)


class PSD_AZ(geometry.SegmentedPath):
    def __init__(self, pointli, profile=None):
        geometry.SegmentedPath.__init__(self, pointli)
        self.profile = profile

    def determineStuff(self):
        self.__adjust()
        self.synm = self.getSynm()
        self.outline = self.__getOutline()
        self.psd_az_loc = self.centerPoint


    def __adjust(self):
        """ Adjust PSD/active zone coordinates so that the PSD/active zone ends 
            exactly on the plasma membrane. 
        """
        # Partition psd/active zone into paths defined by the intersections with
        # the plasma membrane, beginning and ending with the projections of the 
        # end nodes of psd/active zone. 
        # Each path will then be completely on one side of plasma membrane.
        pathli = [] 
        pathli.append(geometry.SegmentedPath())
        pathli[0].append(self[0].projectOnClosedPath(self.profile.path)[0])
        for n in range(0, len(self)-1): 
            pathli[-1].append(self[n])
            for d in range(len(self.profile.path)-1):
                p = geometry.segmentIntersection(self[n], self[n+1],
                                self.profile.path[d], self.profile.path[d+1])
                if p:
                    break
            if p: 
                pathli[-1].append(p)   
                pathli.append(geometry.SegmentedPath())
                pathli[-1].append(p)     
        pathli[-1].append(self[-1])
        pathli[-1].append(self[-1].projectOnClosedPath(self.profile.path)[0])
        # Now, look for the longest path. This is assumed to be the intended
        # psd. (Perhaps area is more relevant. However, this is messier because
        # we need to determine the part of dm enclosed by path for each path.)
        maxLength = 0
        for path in pathli:
            L = path.length
            if L > maxLength:
                longestPath = path
                maxLength = L
        del self[:]
        self.extend(longestPath)
        # Now, orient the PSD so that it runs in the same direction as the
        # profile border
        self.orientToPath(self.profile.path)
        return self

    def getSynm(self):
        """ Return synaptic membrane (i e the part of the plasma membrane
            bounded by the PSD/active zone). Assume that the PSD/active zone is 
            adjusted so as to end exactly on plasma membrane, and
            that the PSD/active zone is oriented in the same direction as the
            plasma membrane (ie node1 <= node2).
        """
        synm = geometry.SegmentedPath()
        p1, node1 = self[0].projectOnClosedPath(self.profile.path)
        p2, node2 = self[-1].projectOnClosedPath(self.profile.path)
        if not None in (node1, node2):
            if node1 == node2: # because of iteratePartial(), we need to 
                               # specifically deal with this special case
                synm.extend([p1, p2])
            else:    
                synm.extend([p for p in self.profile.path.iteratePartial(
                                                                     node1+1,
                                                                     node2+1)])
                synm.insert(0, p1)
                synm.append(p2)
                if synm.length > (self.profile.path.perimeter / 2):
                    synm = geometry.SegmentedPath()
                    synm.extend([p for p in 
                                         self.profile.path.iteratePartial(
                                                                      node2+1,
                                                                      node1+1)])
                    synm.insert(0, p2)
                    synm.append(p1)
                    synm.reverse()
        if len(synm) == 0 or synm.length == 0:
            raise ProfileError(self.profile,
                                "Could not determine synaptic membrane")
        # It appears that sometimes the first or last nodes get duplicated, 
        # resulting in a zero vector and division by zero when projecting on 
        # synm. Simply checking for the duplicated node and removing it if 
        # present obviously works, but perhaps I should delve deeper into the 
        # problem. The node should not be duplicated in the first place. 
        # - 12/15/2003
        try:
            if synm[0] == synm[1]:
                del synm[0]
        except IndexError:
            pass        
        try:
            if synm[-1] == synm[-2]:
                del synm[-1]
        except IndexError:
            pass
        return synm


    def __getOutline(self):
        """ Make a closed-path PSD by concatenating the PSD and the PoSM.
            Assume that PSD and synm share end nodes and that they are 
            oriented in the same direction. 
        """
        pol = geometry.SegmentedPath()
        pol.extend(self.synm)
        pol.reverse()
        pol.extend(self[1:-1])
        return pol 
                          
# End of class PSD_AZ

class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.warnflag = False
        self.opt = opt
        self.errflag = False             

    def process(self, opt):
        """ Parse profile data from a file and determine distances
        """
        try:
            self.parse(opt)
            self.checkPaths()
            sys.stdout.write("Processing profile...\n")
            self.spatial_resolution_in_pixels = geometry.toPixelUnits(
                                                    self.opt.spatial_resolution,
                                                    self.pixelwidth)                             
            self.shell_width_in_pixels = geometry.toPixelUnits(
                                                    self.opt.shell_width,
                                                    self.pixelwidth)
            for psd in self.psd_az_li:
                psd.determineStuff()
            self.totalSynm = self.getTotalSynm()
            self.posloc = self.psd_az_li[0].psd_az_loc
            self.pli.determineStuff()
            if self.opt.useGrid:
                self.gli.determineStuff()
            if self.opt.useRandom:
                self.randomli.determineStuff()
            self.getInterDistlis()
            self.getClusters()
            self.getMonteCarlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
            if self.opt.outputs['individual profiles']:
                self.saveResults()
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True
            
    def getMonteCarlo(self):
        if self.opt.run_monte_carlo == False:
            self.mcli = []
        else:
            sys.stdout.write("Running Monte Carlo simulations...\n")
            self.mcli = self.runMonteCarlo()

    def getInterDistlis(self):
        if not self.opt.determine_interpoint_dists:
            return
        if not True in [val for key, val in
                        self.opt.interpoint_relations.items()
                                if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if self.opt.interpoint_relations["particle - particle"]:
            self.pp_distli, self.pp_latdistli = self.getSameInterpointDistances(
                                                                    self.pli)
        if (self.opt.useRandom and
            self.opt.interpoint_relations["random - particle"]):
            self.rp_distli, self.rp_latdistli = \
                self.getInterpointDistances2(self.randomli, self.pli)


    def getClusters(self):
        if not (self.opt.determine_clusters and
                self.opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.clusterli = self.determineClusters(self.pli)
        self.processClusters(self.clusterli)

    def getSameInterpointDistances(self, pointli):
        dli = []
        latdli = []
        for i in range(0, len(pointli)):
            if self.opt.stop_requested:
                return [], []
            if self.opt.interpoint_dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if self.opt.interpoint_shortest_dist:
                        dli.append(pointli[i].dist(pointli[j]))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(pointli[i].lateralDistToPoint(
                                                                   pointli[j],
                                                                   self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].determineNearestNeighbour(pointli))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].determineNearestLateralNeighbour(
                                                            pointli))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def getInterpointDistances2(self, pointli, pointli2=[]):
        dli = []
        latdli = []
        for i, p in enumerate(pointli):
            if self.opt.stop_requested:
                return [], []
            if self.opt.interpoint_dist_mode == 'all':
                for p2 in pointli2:
                    if self.opt.interpoint_shortest_dist:
                        dli.append(p.dist(p2))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(p.lateralDistToPoint(p2, self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(p.determineNearestNeighbour(pointli2))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(p.determineNearestLateralNeighbour(pointli2))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def runMonteCarlo(self):

        def isValid(p):
            def isWithinOrAssociatedWithProfile():
                # Better to use this as p.distToPath has already been computed
                if ((p.distToPath >= 0 or
                    (self.opt.monte_carlo_strict_location == False and
                     abs(p.distToPath) <= self.spatial_resolution_in_pixels))):
                    return True
                else:
                    return False

            if p in mcli[n]["pli"] or p.isWithinHole:
                return False
            p.distToPath = p.getDistToPath()
            if p.distToPath < 0 and abs(p.distToPath) > border:
                return False
            if self.opt.monte_carlo_simulation_window == "profile + shell":
                return True    # Points outside shell have been discarded above
            if (self.opt.monte_carlo_simulation_window == "profile" and
                isWithinOrAssociatedWithProfile()):
                return True
            if (self.opt.monte_carlo_simulation_window ==
                "profile - synaptic zone" and
                isWithinOrAssociatedWithProfile() and not 
                p.getIsWithinPSD_AZ):  # only discard points strictly
                return True                   # within PSD/AZ regardless of
                                              # self.opt.monte_carlo_strict_location
            if (self.opt.monte_carlo_simulation_window ==
                "profile + shell - synaptic zone" and not
                p.getIsWithinPSD_AZ):
                return True
            if (self.opt.monte_carlo_simulation_window == "synaptic zone"):
                if self.opt.monte_carlo_strict_location == True:
                    if p.getIsWithinPSD_AZ:
                        return True
                elif p.getIsAssociatedWithPSD_AZ:
                    return True
            return False

        sys.stdout.write("-> Generating simulated points:\n   ")
        pli = self.pli
        box = self.path.boundingBox()
        border = self.shell_width_in_pixels
        if "profile + shell" in self.opt.monte_carlo_simulation_window:
            numpoints = len(pli)
        elif "profile" in self.opt.monte_carlo_simulation_window:
            if self.opt.monte_carlo_strict_location:
                numpoints = len([p for p in pli if p.isWithinProfile])
            else:
                numpoints = len([p for p in pli if p.isAssociatedWithProfile])
        elif self.opt.monte_carlo_simulation_window == "synaptic zone":
            if self.opt.monte_carlo_strict_location:
                numpoints = len([p for p in pli if p.getIsWithinPSD_AZ])
            else:
                numpoints = len([p for p in pli if p.getIsAssociatedWithPSD_AZ])
        else:
            return []
        if "- synaptic zone" in self.opt.monte_carlo_simulation_window:
            numpoints -= len([p for p in pli if p.getIsWithinPSD_AZ])
        mcli = []
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []
            dotProgress(n)
            mcli.append({"pli": [],
                         "simulated - simulated": {"dist": [], "latdist": []},
                         "simulated - particle": {"dist": [], "latdist": []},
                         "particle - simulated": {"dist": [], "latdist": []},
                         "clusterli": []})
            for i in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - border),
                                       int(box[1].x + border) + 1)
                    y = random.randint(int(box[0].y - border),
                                       int(box[2].y + border) + 1)
                    p = Point(x, y, ptype="simulated", profile=self)
                    if isValid(p):
                        break
                # escape the while loop when a valid simulated point is found
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determineStuff()
        sys.stdout.write("\n-> Determining simulated interpoint distances:\n   ")
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []
            dotProgress(n)
            if self.opt.interpoint_relations["simulated - simulated"]:
                    distlis = self.getSameInterpointDistances(mcli[n]["pli"])
                    mcli[n]["simulated - simulated"]["dist"].append(distlis[0])
                    mcli[n]["simulated - simulated"]["latdist"].append(distlis[1])
            if self.opt.interpoint_relations["simulated - particle"]:
                distlis = self.getInterpointDistances2(mcli[n]["pli"], pli)
                mcli[n]["simulated - particle"]["dist"].append(distlis[0])
                mcli[n]["simulated - particle"]["latdist"].append(distlis[1])
            if self.opt.interpoint_relations["particle - simulated"]:
                distlis = self.getInterpointDistances2(pli, mcli[n]["pli"])
                mcli[n]["particle - simulated"]["dist"].append(distlis[0])
                mcli[n]["particle - simulated"]["latdist"].append(distlis[1])
        sys.stdout.write("\n-> Determining simulated clusters:\n   ")
        if self.opt.determine_clusters:
            for n, li in enumerate(mcli):
                if self.opt.stop_requested:
                    return []
                dotProgress(n)
                mcli[n]["clusterli"] = self.determineClusters(li["pli"])
                self.processClusters(mcli[n]["clusterli"])
        sys.stdout.write("\n")
        return mcli


    def processClusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.convexHull = geometry.convexHullGraham(c)
            c.distToPath = c.convexHull.centroid.perpendDistClosedPath(self.path)
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.nearestCluster = ClusterData()
            if len(clusterli) == 1:
                c.distToNearestCluster = -1
                return
            c.distToNearestCluster = sys.maxint
            for c2 in clusterli:
                if c2 is not c:
                    d = c.lateralDist(c2, self.path)
                    if  d < c.distToNearestCluster:
                        c.distToNearestCluster = d
                        c.nearestCluster = c2


    def determineClusters(self, pointli):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than self.opt.within_cluster_dist from at least one
            other point in the cluster
        """
        clusterli = []
        for p1 in pointli:
            if self.opt.stop_requested:
                return []
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.toPixelUnits(
                                                    self.opt.within_cluster_dist,
                                                    self.pixelwidth):
                    if p2.cluster != None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        return clusterli


    def parse(self, opt):
        """ Parse profile data from input file 
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = readFile(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        self.psd_az_li = []
        self.holeli = []
        while li:
            s = li.pop(0).replace("\n","").strip()
            if s.split(" ")[0].upper() == "IMAGE":
                self.src_img = s.split(" ")[1]
            elif s.split(" ")[0].upper() == "PROFILE_ID":
                try:
                    self.ID = s.split(" ")[1]
                except IndexError, ValueError:
                    ProfileWarning(self, "Profile ID not defined or invalid")
            elif s.split(" ")[0].upper() == "COMMENT":
                try:
                    self.comment = s.split(" ", 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(" ")[0].upper() == "PIXELWIDTH":
                try: 
                    self.pixelwidth = float(s.split(" ")[1])
                    self.metric_unit = s.split(" ")[2]
                except ValueError:
                    raise ProfileError(self, 
                                       "PIXELWIDTH is not a valid number")
            elif s.upper() in ("PROFILE_BORDER", "PLASMA_MEMBRANE"):
                self.path = PlasmaMembraneData(self.getCoords(li, "path"),
                                                profile=self)
            elif s.upper() == "PSD_OR_ACTIVE_ZONE":
                self.psd_az_li.append(PSD_AZ(self.getCoords(li, "PSD_AZ"),
                                              profile=self))
            elif s.upper() == "PROFILE_HOLE":
                self.path.addHole(geometry.SegmentedPath(self.getCoords(li, "hole")))
            elif s.upper() == "POINTS":
                #pli = self.getCoords(li, "point")
                self.pli = PointList(self.getCoords(li, "point"),
                                     "point", profile=self)
                #for p in pli:
                #    self.pli.append(Point(p.x, p.y, ptype="point",
                #                    profile=self))
            elif s.upper() == "GRID":
                self.gridli = PointList(self.getCoords(li, "grid"),
                                           "grid", profile=self)
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(self.getCoords(li, "random"),
                                             "random", profile=self)
            elif s[0] != "#":          # unless specifically commented out           
                ProfileWarning(self, "Unrecognized string '" + s + 
                                     "' in input file")
        # Now, let's see if everything was found
        self.checkParsedData(opt)

    def checkParsedData(self, opt):
        """ See if the synapse data was parsed correctly, and print info on the
            parsed data to standard output.            
        """
        self.checkVarDefault(self, 'src_img', "Source image", "N/A")
        self.checkVarDefault(self, 'ID', "Profile ID", "N/A")
        self.checkVarDefault(self, 'comment', "Comment", "")
        self.checkVarVal(self, 'metric_unit', "Metric unit", 'metric_unit', opt)
        self.checkRequiredVar(self, 'pixelwidth', "Pixel width", self.metric_unit)
        self.checkVarDefault(self, 'postsynProfile', "Postsynaptic profile", "N/A")
        self.checkVarDefault(self, 'presynProfile', "Presynaptic profile", "N/A")
        self.checkListVar(self, 'path', 'Plasma membrane', 'nodes', 3)
        self.checkTableVar(self, 'psd_az_li', "Postsynaptic density/active zone",
                             "Postsynaptic densities/active zones", 1, 2)
        self.checkListVar(self, 'pli', 'Particles', '', 0)
        self.checkTableVar(self.path, 'holeli', "Hole", "Holes", 0, 2)
        self.checkVarExists(self, 'gridli', "Grid", 'useGrid', opt)
        self.checkVarExists(self, 'randomli', "Random points", 'useRandom', opt)
        if not self.path.isSimplePolygon():
            raise ProfileError(self, "Profile is not a simple polygon")
        for n, h in enumerate(self.path.holeli):
            if not h.isSimplePolygon():
                raise ProfileError(self, 
                                   "Profile hole %d is not a simple polygon" 
                                    % (n+1))
            if not h.isWithinPolygon(self.path):
                raise ProfileError(self, 
                                   "Profile hole %d is not (completely) "
                                   "within profile" % (n+1))                
            for n2, h2 in enumerate(self.path.holeli[n+1:]):
                if h.overlapsPolygon(h2):
                    raise ProfileError(self, 
                                       "Profile hole %d overlaps with hole %d "
                                       % (n+1, n+n2+2))                                    
                
    def checkRequiredVar(self, parent, var_to_check, var_str, post_str):
        """ Confirm that a required variable exists; else, raise ProfileError.
        """
        if not hasattr(parent, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n"
                             % (var_str, parent.__dict__[var_to_check], post_str))

    def checkListLen(self, var, min_len):
        """ Returns True if var is a list and has at least min_len elements,
            else False
        """
        return isinstance(var, list) and len(var) >= min_len

    def checkListVar(self, parent, var_to_check, var_str, post_str, min_len):
        """ Confirms that parent has a var_to_check that is a list and has at
            least min_len elements; if var_to_check does not exist and
            min_len <= 0, assigns an empty list to var_to_check. Else, raise a
            ProfileError.
        """
        if not hasattr(parent, var_to_check):
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file"
                                         % var_to_check)
            else:
                parent.__dict__[var_to_check] = []
        elif not self.checkListLen(parent.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n"
                         % (var_str, post_str,
                            len(parent.__dict__[var_to_check])))


    def checkTableVar(self, parent, var_to_check, var_str_singular, var_str_plural,
                        min_len_1, min_len_2):
        """ Confirms that var_to_check exists, is a list and has at least
            min_len_1 elements, and that each of these has at least min_len_2
            subelements; if var_to_check does not exist and min_len_1 <= 0,
            assigns an empty list to var_to_check. Else, raise ProfileError.
        """
        if not hasattr(parent, var_to_check):
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file"
                                         % var_str_plural)
            else:
                parent.__dict__[var_to_check] = []
        elif not self.checkListLen(parent.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file"
                               % var_str_plural.lower())
        else:
            for element in parent.__dict__[var_to_check]:
                if not self.checkListLen(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates"
                                       % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural,
                                         len(parent.__dict__[var_to_check])))

    def checkVarExists(self, parent, var_to_check, var_str, optflag, opt):
        """ Checks for consistency between profiles with respect to the
            existence of var_to_check (i.e., var_to_check must be present
            either in all profiles or in none).

            If optflag is not set (i.e., this is the first profile), then
            set optflag to True or False depending on the existence of
            var_to_check. If optflag is already set (for consequent profiles),
            var_to_check must (if optflag is True) or must not (if optflag is
            False) exist. If not so, raise ProfileError.
        """
        if not hasattr(opt, optflag):
            if hasattr(self, var_to_check):
                opt.__dict__[optflag] = True
            else:
                opt.__dict__[optflag] = False
        if opt.__dict__[optflag]:
            if hasattr(parent, var_to_check):
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s not found in input file" % var_str)
        elif hasattr(parent, var_to_check):
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def checkVarVal(self, parent, var_to_check, var_str, optvar, opt):
        """ Checks for consistency between profiles with respect to the
            value of var_to_check (i.e., var_to_check must be present and
            have equal value in all profiles).

            If optvar is not set (i.e., this is the first profile), then
            set optflag to the value of var_to_check. If optvar is already set
            (for consequent profiles), the value of var_to_check must be equal
            to that of optvar. If not so, raise ProfileError.
        """
        if not hasattr(parent, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(opt, optvar):
            opt.__dict__[optvar] = parent.__dict__[var_to_check]
        elif parent.__dict__[var_to_check] == opt.__dict__[optvar]:
            pass # really no point in pointing out that it's ok
            #sys.stdout.write("  %s: %s\n"
            #                 % (var_str, parent.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, parent.__dict__[var_to_check],
                                  opt.__dict__[optvar]))

    def checkVarDefault(self, parent, var_to_check, var_str, default=""):
        """ Checks if var_to_check exists; if not, assign the default value
        to var_to_check. Never raises a ProfileError.
        """
        if not hasattr(parent, var_to_check):
            parent.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str,
                                         parent.__dict__[var_to_check]))
    def checkPaths(self):
        """ Make sure that plasma membrane, PSDs/active zones and holes 
            do not intersect with themselves
        """
        def checkPath(path, s):
            for n1 in range(0, len(path)-3):
                for n2 in range(0, len(path)-1):
                    if n1 not in (n2, n2+1) and n1+1 not in (n2, n2+1):
                        if geometry.segmentIntersection(path[n1], path[n1+1],
                                               path[n2], path[n2+1]):
                            raise ProfileError(self, "%s invalid (crosses itself)" % s)
            return True
                    
        while self.path[-1] == self.path[0]:
            del self.path[-1]
        checkPath(self.path, "Plasma membrane")
        for path in self.psd_az_li:
            checkPath(path, "PSD/active zone")
        for path in self.path.holeli:
            while path[-1] == path[0]:
                del path[-1]
            checkPath(path, "Hole")
        sys.stdout.write("  Paths are ok.\n")
 
 
    def getTotalSynm(self):
        """ Construct a path comprising the plasma membrane of all 
            PSDs/active zones and perforations.
            Assume well-behaved PSDs/active zones.
        """

        def distToLeftEndnode(p):
            path = geometry.SegmentedPath()
            project, seg_project = p.projectOnPathOrEndnode(self.path)
            path.extend([self.path[0], project])
            for n in range(1, seg_project):
                path.insert(len(path)-1, self.path[n])
            return path.length

        left_mindist = right_mindist = self.path.length
        for p in itertools.chain(*[[psd_az_[0], psd_az_[-1]]
                                    for psd_az_ in self.psd_az_li]):
            left_d = distToLeftEndnode(p)
            right_d = self.path.length - left_d
            if left_d <= left_mindist:
                left_mindist = left_d
                left_p = p
            elif right_d <= right_mindist:
                right_mindist = right_d
                right_p = p
        pseudoPSD_AZ = PSD_AZ(geometry.SegmentedPath([left_p, right_p]), self)
        return pseudoPSD_AZ.getSynm()
        
        

    def getCoords(self, strli, coordType=""):
        """ Pop point coordinates from list strli.
            When an element of strli is not a valid point,
            a warning is issued.
        """
        pointli = []
        s = strli.pop(0).replace("\n","").replace(" ","").strip()
        while s != "END":
            try:
                p = geometry.Point(float(s.split(",")[0]), float(s.split(",")[1]))
                if pointli and (p == pointli[-1] or 
                                (coordType in ('point', 'random')
                                 and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coordType, p))
                else:
                    pointli.append(Point(p.x, p.y, ptype=coordType,
                                   profile=self))
            except ValueError:
                if s[0] != "#":
                    ProfileWarning(self, "'%s' not valid %s coordinates" 
                                   % (s, coordType))
                else:
                    pass 
            s = strli.pop(0).replace("\n","").strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]): 
            del pointli[-1]                                            
        return pointli        

    def saveResults(self):
        """ Output results from a single profile to file
        """ 
        
        def m(x):
            try:
                return geometry.toMetricUnits(x, self.pixelwidth)
            except ZeroDivisionError:
                return None
                
        def m2(x): 
            try:
                return geometry.toMetricUnits(x, self.pixelwidth**2) # for area units...
            except ZeroDivisionError:
                return None
       
        def fwrite(*args):
            f.writerow(args)
                    
        
        try:
            self.outputfn = os.path.join(self.opt.output_dir,
                                         os.path.basename(self.inputfn)
                                         + self.opt.output_filename_suffix
                                         + self.opt.output_filename_ext)

            if (os.path.exists(self.outputfn) and
                self.opt.action_if_output_file_exists == 'enumerate'):
                    self.outputfn = enumFilename(self.outputfn, 2)
            sys.stdout.write("Writing to '%s'...\n" % self.outputfn)
            if self.opt.output_file_format == "csv":
                csv_format = { 'dialect' : 'excel', 'lineterminator' : '\n'}
                if self.opt.csv_delimiter == 'tab':
                    csv_format['delimiter'] = '\t'
                f = unicode_csv.writer(file(self.outputfn, "w"),
                                        **self.opt.csv_format)
            elif self.opt.output_file_format == 'excel':
                import xls
                f = xls.writer(self.outputfn)
            fwrite("Table 1. Profile-centric data")
            fwrite("Source image:", self.src_img)
            fwrite("Profile ID:", self.ID)
            fwrite("Comment:", self.comment) 
            fwrite("Pixel width:", tostr(float(self.pixelwidth), 2), 
                                   self.metric_unit)
            fwrite("Number of points (total):", len(self.pli))            
            fwrite("Number of random points (total):", len(self.randomli))
            fwrite("Table 2. Point-centric data")                                   
            columnheadings = ["Point number (as appearing in input file)",
                              "Point coordinates (in pixels)"]
            fwrite(*columnheadings) 
            f.writerows([[n+1,
                          str(p)] 
                          for n, p in enumerate(self.pli)])     
            fwrite("Table 3. Random point-centric data")                                   
            columnheadings = ["Random point number (as appearing in input file)",
                              "Random point coordinates (in pixels)"]
            fwrite(*columnheadings) 
            f.writerows([[n+1,
                          str(r)] 
                          for n, r in enumerate(self.randomli)])                               
            f.close()
        except IOError:
            raise ProfileError(self, "Unable to write to file '%s'" 
                               % self.outputfn)
            return 0
        sys.stdout.write("Done.\n")
        return 1


# end of class Profile        
        
class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 200   # Skip points farther than this from the
                                # postsynaptic element
        self.outputs = {
                        'profile summary': True, 'point summary': True,
                        'random summary': True, 'session summary': True,
                        'individual profiles': False
                        }
        self.output_file_format = "excel"
        self.output_filename_ext = ".xls"
        self.input_filename_ext = ".pds"
        self.output_filename_suffix = ''
        self.output_filename_other_suffix = ''
        self.output_filename_date_suffix = True
        self.output_filename_use_other_suffix = False
        self.csv_delimiter = 'comma'
        self.action_if_output_file_exists = 'overwrite'
        self.output_dir = ''
        self.determine_clusters = False
        self.within_cluster_dist = 50
        self.run_monte_carlo = False
        self.monte_carlo_runs = 99
        self.monte_carlo_simulation_window = "profile + shell"
        self.monte_carlo_strict_location = False
        self.determine_interpoint_dists = False
        self.interpoint_dist_mode = 'nearest neighbour'
        self.interpoint_relations = {
                                    'particle - particle': True,
                                    'random - particle': True,
                                    'particle - simulated': False,
                                    'simulated - particle': False,
                                    'simulated - simulated': False
                                    }
        self.interpoint_shortest_dist = True
        self.interpoint_lateral_dist = False

    def reset(self):
        """ Resets all options to default, and removes those that are not
            set in __init__().
        """
        self.__dict__ = {}
        self.__init__()
# end of class OptionData

class ProfileError(exceptions.Exception):
    def __init__(self, profile, msg):
        self.args = (profile, msg + ".")

def ProfileWarning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True      

def ProfileMessage(profile, msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)
