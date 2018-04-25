import sys
import random
import itertools
from . import geometry
from . import file_io


# Convenience functions

def dot_progress(line_length=80, char='.', reset=False):
    """Simple progress indicator on sys.stdout"""
    if not hasattr(dot_progress, 'counter'):
        dot_progress.counter = 0
    if reset:
        dot_progress.counter = 0
        sys.stdout.write('\n')
    dot_progress.counter += 1
    sys.stdout.write(char)
    if dot_progress.counter == line_length:
        dot_progress.counter = 0
        sys.stdout.write('\n')


def lazy_property(fn):
    """Decorator that makes a property lazily evaluated.
       From https://stevenloria.com/lazy-properties/.
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


#
# Classes
#


class Point(geometry.Point):
    def __init__(self, x=None, y=None, ptype='', profile=None):
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile is not None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.discard = False
        self.ptype = ptype
        self.cluster = None
        self.is_within_psd = None
        self.is_associated_with_psd = None
        self.associated_psd = None
        self.lateral_dist_psd = None
        self.norm_lateral_dist_psd = None
        self.perpend_dist_psd = None
        self.nearest_psd = None
        self.nearest_neighbour_dist = None
        self.nearest_neighbour_point = geometry.Point()
        self.nearest_lateral_neighbour_dist = None
        self.nearest_lateral_neighbour_point = geometry.Point()
        self.nearest_neighbour = geometry.Point()

    def determine_stuff(self):
        """Determine general stuff for a point, including distance to path.
         Also mark the point for discarding if it is not valid.
        """

        def mark_to_discard(msg):
            if self.ptype == 'particle':  # don't warn if random points
                profile_message("Discarding particle at %s: %s" % (self, msg))
            self.discard = True
            self.profile.n_discarded[self.ptype] += 1
            return

        if self.is_within_hole:
            mark_to_discard("Located within a profile hole")
            return
        if not (self.is_within_profile or self.is_within_shell):
            mark_to_discard("Located outside the shell")
            return
        # This is to force the computation of this lazy property here
        __ = self.is_associated_with_path
        self.perpendicular_dist_syn(self.profile.path, self.profile.total_synm)
        self.get_psd_association()
        self.get_lateral_dist_psd()

    @lazy_property
    def dist_to_path(self):
        """Return distance to profile border"""
        _dist_to_path = self.perpend_dist_closed_path(self.profile.path)
        if not self.is_within_profile:
            _dist_to_path = -_dist_to_path
        return _dist_to_path

    @lazy_property
    def is_within_hole(self):
        """Determine whether self is inside a profile hole"""
        within_hole = False
        for h in self.profile.holeli:
            if self.is_within_polygon(h):
                within_hole = True
            else:
                within_hole = False
        return within_hole

    @lazy_property
    def is_within_profile(self):
        """Determine whether self is inside profile, excluding holes"""
        if self.is_within_polygon(self.profile.path) and not self.is_within_hole:
            return True
        else:
            return False

    @lazy_property
    def is_within_shell(self):
        """Determine whether self is within shell"""
        return (not self.is_within_profile and
                abs(self.dist_to_path) < geometry.to_pixel_units(
                    self.opt.shell_width, self.profile.pixelwidth))

    @lazy_property
    def is_associated_with_path(self):
        """Determine whether self is associated with the profile
        border, i e, is within a distance of it that is less than
        the spatial resolution"""
        if (abs(self.dist_to_path) <= geometry.to_pixel_units(
                self.opt.spatial_resolution, self.profile.pixelwidth)):
            return True
        else:
            return False

    @lazy_property
    def is_associated_with_profile(self):
        """Determine whether self is within the profile or
        associated with the profile border"""
        if self.is_within_profile or self.is_associated_with_path:
            return True
        else:
            return False

    def get_lateral_dist_psd(self):
        """Determine nearest PSD and metric and normalized
        lateral distance to this PSD.
        """
        mindist = sys.maxsize
        nearest_psd = None
        for psd in self.profile.psd_li:
            d = (self.lateral_dist_syn(self.profile.path, psd.synm) - (psd.synm.length() / 2))
            if d < mindist:
                mindist = d
                nearest_psd = psd
        if nearest_psd is None:
            raise ProfileError(self.profile,
                               "could not determine lateral distance to a "
                               "PSD of point at %s" % self)
        mindist = self.lateral_dist_syn(self.profile.path, nearest_psd.synm)
        normdist = mindist / (nearest_psd.synm.length() / 2)
        self.lateral_dist_psd = mindist
        self.norm_lateral_dist_psd = normdist
        self.nearest_psd = nearest_psd

    def lateral_dist_syn(self, path, sm):
        """Determine lateral distance to center of synaptic membrane
        sm, which is a subpath of path.
        """
        return self.lateral_dist_to_point(sm.center_point(), path)

    def perpendicular_dist_syn(self, pm, sm):
        """Determine perpendicular distance of self to synaptic
        membrane sm; if the projection is not on sm, return distance
        to the nearest endnode of sm.
        """
        d = self.perpend_dist(sm)
        if d is None:
            d = min(abs(self.dist(sm[0])), abs(self.dist(sm[-1])))
        if self.is_within_polygon(pm):
            self.perpend_dist_psd = d
        else:
            self.perpend_dist_psd = -d

    def get_psd_association(self):
        """ Determine whether self is within or associated with a PSD,
        and the associated PSD, if any.
        """
        is_within_psd = False
        is_associated_with_psd = False
        associated_psd = None
        mindist = sys.maxsize
        for psd in self.profile.psd_li:
            if self.is_within_polygon(psd.outline):
                is_within_psd = True
                is_associated_with_psd = True
                associated_psd = psd
                break
            dist = self.perpend_dist_closed_path(psd.outline, dont_care_if_on_or_off_seg=True)
            if dist <= geometry.to_pixel_units(self.profile.opt.spatial_resolution,
                                               self.profile.pixelwidth):
                is_associated_with_psd = True
                if dist < mindist:
                    associated_psd = psd
                    mindist = dist
        self.is_within_psd = is_within_psd
        self.is_associated_with_psd = is_associated_with_psd
        self.associated_psd = associated_psd

    def get_nearest_neighbour(self, pointli):
        """Determine distance to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        minp = Point()
        for p in pointli:
            # I will instead exclude non-desired points from the
            # supplied point list *before* calling this function
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxsize):
            return None
        else:
            self.nearest_neighbour_dist = mindist
            self.nearest_neighbour_point = minp
            return self.nearest_neighbour_dist

    def get_nearest_lateral_neighbour(self, pointli):
        """Determine distance along profile border to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        minp = Point()
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, self.profile.path)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxsize):
            return None
        else:
            self.nearest_lateral_neighbour_dist = mindist
            self.nearest_lateral_neighbour_point = minp
            return self.nearest_lateral_neighbour_dist


class ProfileBorderData(geometry.SegmentedPath):
    def __init__(self, pointlist=None, profile=None):
        if pointlist is None:
            pointlist = []
        geometry.SegmentedPath.__init__(self, pointlist)
        self.profile = profile


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        super().__init__()
        try:
            self.extend([Point(p.x, p.y, ptype, profile) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError("not a list of Point elements")


class ClusterData(list):
    def __init__(self, pointli=None):
        super().__init__()
        if pointli is None:
            pointli = []
        try:
            self.extend([Point(p.x, p.y) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError("not a point list")
        self.convex_hull = geometry.SegmentedPath()

    def lateral_dist_to_cluster(self, c2, border):
        """Determine lateral distance to a cluster c2 along profile
        border.
        """
        path = geometry.SegmentedPath()
        c2_project, c2_seg_project = c2.convex_hull.centroid().\
            project_on_closed_path(border)
        project, seg_project = self.convex_hull.centroid().\
            project_on_closed_path(border)
        path.extend([project, c2_project])
        if c2_seg_project < seg_project:
            path.reverse()
        for n in range(min(c2_seg_project, seg_project) + 1,
                       max(c2_seg_project, seg_project)):
            path.insert(len(path) - 1, border[n])
        length = path.length()
        return min(length, border.perimeter() - length)


class Psd(geometry.SegmentedPath):
    def __init__(self, pointli, profile):
        geometry.SegmentedPath.__init__(self, pointli)
        self.profile = profile
        self.synm = geometry.SegmentedPath()
        self.outline = geometry.SegmentedPath()
        self.psd_loc = geometry.Point()

    def determine_stuff(self):
        self.__adjust()
        self.synm = self.get_synm()
        self.outline = self.__get_outline()
        self.psd_loc = self.center_point()

    def __adjust(self):
        """ Adjust PSD coordinates so that the PSD ends exactly
         on the plasma membrane.
        """
        # TODO: Revise to deal better with "bad" PSD outlines
        # Partition psd into paths defined by the intersections with
        # the plasma membrane, beginning and ending with the projections
        # of the end nodes of psd. Each path will then be completely on
        # one side of plasma membrane.
        pathli = [geometry.SegmentedPath()]
        pathli[0].append(self[0].project_on_closed_path(self.profile.path)[0])
        p = Point()
        for n in range(0, len(self)-1): 
            pathli[-1].append(self[n])
            for d in range(len(self.profile.path)-1):
                p = geometry.segment_intersection(self[n],
                                                  self[n+1],
                                                  self.profile.path[d],
                                                  self.profile.path[d+1])
                if p:
                    break
            if p: 
                pathli[-1].append(p)   
                pathli.append(geometry.SegmentedPath())
                pathli[-1].append(p)     
        pathli[-1].append(self[-1])
        pathli[-1].append(self[-1].project_on_closed_path(self.profile.path)[0])
        # Now, look for the longest path. This is assumed to be the intended
        # psd. (Perhaps area is more relevant. However, this is messier because
        # we need to determine the part of dm enclosed by path for each path.)
        # Update 2016-10-14: the length of the path is here defined as its
        # number of vertices.
        max_length = 0
        old_max_length = 0
        longest_path = geometry.SegmentedPath()
        for path in pathli:
            if len(path) >= max_length:
                longest_path = path
                old_max_length = max_length
                max_length = len(path)
        if max_length < 3 or max_length == old_max_length:
            raise ProfileError(self.profile, "Could not determine postsynaptic density")
        del self[:]
        self.extend(longest_path)
        # Now, orient the PSD so that it runs in the same direction as the
        # profile border
        self.orient_to_path(self.profile.path)
        return self

    def get_synm(self):
        """ Return synaptic membrane (i e the part of the plasma
        membrane bounded by the PSD). Assume that the PSD is
        adjusted so as to end exactly on plasma membrane membrane,
        and that the PSD is oriented in the same direction
        as the plasma membrane (ie node1 <= node2).
        """
        synm = geometry.SegmentedPath()
        p1, node1 = self[0].project_on_closed_path(self.profile.path)
        p2, node2 = self[-1].project_on_closed_path(self.profile.path)
        if None not in (node1, node2):
            if node1 == node2:  # because of iterate_partial(), we need to
                                # specifically deal with this special case
                synm.extend([p1, p2])
            else:    
                synm.extend([p for p in self.profile.path.iterate_partial(node1+1, node2+1)])
                synm.insert(0, p1)
                synm.append(p2)
                if synm.length() > (self.profile.path.perimeter() / 2):
                    synm = geometry.SegmentedPath()
                    synm.extend([p for p in self.profile.path.iterate_partial(node2+1, node1+1)])
                    synm.insert(0, p2)
                    synm.append(p1)
                    synm.reverse()
        if len(synm) == 0 or synm.length == 0:
            raise ProfileError(self.profile, "Could not determine synaptic membrane")
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

    def __get_outline(self):
        """ Make a closed-path PSD by concatenating the PSD and the synm.
            Assume that PSD and synm share end nodes and that they are 
            oriented in the same direction. 
        """
        pol = geometry.SegmentedPath()
        pol.extend(self.synm)
        pol.reverse()
        pol.extend(self[1:-1])
        return pol
# End of class PsdAz


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.src_img = None
        self.opt = opt
        self.holeli = []
        self.psd_li = []
        self.pli = []
        self.randomli = []
        self.mcli = []
        self.clusterli = []
        self.pp_distli, self.pp_latdistli = [], []
        self.rp_distli, self.rp_latdistli = [], []
        self.n_discarded = {'particle': 0, 'random': 0}
        self.comment = ''
        self.pixelwidth = None
        self.metric_unit = ''
        self.posloc = geometry.Point()
        self.negloc = geometry.Point()
        self.total_synm = geometry.SegmentedPath()
        self.perimeter = None
        self.feret = None
        self.warnflag = False
        self.errflag = False

    def process(self):
        """ Parse profile data from a file and determine distances
        """
        try:
            sys.stdout.write("Processing profile...\n")
            self.__parse()
            self.__check_paths()
            for psd in self.psd_li:
                psd.determine_stuff()
            self.total_synm = self.get_total_synm()
            self.posloc = self.psd_li[0].psd_loc
            self.__compute_stuff()
            if self.opt.determine_interpoint_dists:
                sys.stdout.write("Determining interparticle distances...\n")
                self.__determine_interdistlis()
            if self.opt.determine_clusters:
                sys.stdout.write("Determining clusters...\n")
                self.clusterli = self.__determine_clusters(self.pli)
            if self.opt.run_monte_carlo:
                sys.stdout.write("Running Monte Carlo simulations...\n")
                self.__run_monte_carlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
        except ProfileError as err:
            sys.stdout.write("Error: %s\n" % err.msg)
            self.errflag = True

    @lazy_property
    def area(self):
        """Determine area of profile, excluding holes"""
        tot_hole_area = sum([h.area() for h in self.holeli])
        return self.path.area() - tot_hole_area

    def contains(self, p):
        """Determine if a point is inside profile, excluding holes."""
        if not p:
            return None
        return p.is_within_profile(self)

    def __compute_stuff(self):
        __ = self.area  # Force computation here
        self.perimeter = self.path.perimeter()
        self.feret = self.path.feret_diameter()
        for p in self.pli:
            p.determine_stuff()
        self.pli = [p for p in self.pli if not p.discard]
        for p in self.randomli:
            p.determine_stuff()
        self.randomli = [p for p in self.randomli if not p.discard]
        for ptype in ('particle', 'random'):
            if ptype == 'random' and not self.opt.use_random:
                continue
            ptypestr = 'particles' if ptype == 'particle' else ptype + ' points'
            sys.stdout.write("  Number of %s discarded: %d\n"
                             % (ptypestr, self.n_discarded[ptype]))

    def __determine_interdistlis(self):
        if True not in [val for key, val in self.opt.interpoint_relations.items()
                        if 'simulated' not in key]:
            return
        if self.opt.interpoint_relations['particle - particle']:
            self.pp_distli, self.pp_latdistli = self.__get_same_interpoint_distances(self.pli)
        if self.opt.use_random and self.opt.interpoint_relations['random - particle']:
            self.rp_distli, self.rp_latdistli = self.__get_interpoint_distances2(self.randomli,
                                                                                 self.pli)

    def __get_same_interpoint_distances(self, pointli):
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
                        latdli.append(pointli[i].lateral_dist_to_point(
                            pointli[j], self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].get_nearest_neighbour(pointli))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].get_nearest_lateral_neighbour(
                        pointli))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def __get_interpoint_distances2(self, pointli, pointli2=None):
        if pointli2 is None:
            pointli2 = []
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
                        latdli.append(p.lateral_dist_to_point(p2, self.path))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(p.get_nearest_neighbour(pointli2))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(p.get_nearest_lateral_neighbour(pointli2))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

# TODO: check simulation windows etc!
    def __run_monte_carlo(self):

        def is_valid(p_candidate):
            def is_within_or_associated_with_profile():
                """Depending on self.opt.monte_carlo_strict_location,
                returns True if within or associated with profile.
                """
                if p_candidate.is_within_profile:
                    return True
                elif self.opt.monte_carlo_strict_location:
                    return False
                # as shell width may be smaller than spatial resolution,
                # we must make sure that the point is also within shell
                elif p_candidate.is_associated_with_profile and p_candidate.is_within_shell:
                    return True
                return False

            # todo: place pt.get_psd_association at some convenient place
            if p_candidate.is_within_hole:
                return False
            if self.opt.monte_carlo_simulation_window == 'profile':
                if is_within_or_associated_with_profile():
                    return True
            elif self.opt.monte_carlo_simulation_window == 'profile + shell':
                if is_within_or_associated_with_profile() or p_candidate.is_within_shell:
                    return True
            elif self.opt.monte_carlo_simulation_window == 'profile - postsynaptic density':
                # only discard points strictly within PSD/AZ regardless of
                # self.opt.monte_carlo_strict_location
                if is_within_or_associated_with_profile() and not p_candidate.is_within_psd:
                    return True
            elif self.opt.monte_carlo_simulation_window == 'profile + shell - postsynaptic density':
                # only discard points strictly within PSD/AZ regardless of
                # self.opt.monte_carlo_strict_location
                if (is_within_or_associated_with_profile() or p_candidate.is_within_shell
                        and not p_candidate.is_within_psd):
                    return True
            elif self.opt.monte_carlo_simulation_window == 'postsynaptic density':
                # print pt.is_within_psd, pt.is_associated_with_psd
                if self.opt.monte_carlo_strict_location and p_candidate.is_within_psd:
                    return True
                elif p_candidate.is_associated_with_psd:
                    return True
            return False

        box = self.path.bounding_box()
        if self.opt.monte_carlo_simulation_window == 'profile' and \
                self.opt.monte_carlo_strict_location:
            border = 0
        else:
            # If shell width is smaller than spatial resolution,
            # the former must be used because all real particles
            # outside the shell have been discarded
            border = geometry.to_pixel_units(min(self.opt.shell_width, self.opt.spatial_resolution),
                                             self.pixelwidth)
        pli = [p for p in self.pli if is_valid(p)]
        numpoints = len(pli)
        mcli = []
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []
            dot_progress(n)
            mcli.append({'pli': [],
                         'simulated - simulated': {'dist': [], 'latdist': []},
                         'simulated - particle': {'dist': [], 'latdist': []},
                         'particle - simulated': {'dist': [], 'latdist': []},
                         'clusterli': []})
            p = Point()
            for __ in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - border),
                                       int(box[1].x + border) + 1)
                    y = random.randint(int(box[0].y - border),
                                       int(box[2].y + border) + 1)
                    p = Point(x, y, ptype='simulated', profile=self)
                    if p not in mcli[n]['pli'] and is_valid(p):
                        break
                # escape the while loop when a valid simulated point is found
                mcli[n]['pli'].append(p)
            for p in mcli[n]['pli']:
                p.determine_stuff()
        if self.opt.determine_interpoint_dists:
            for n in range(0, self.opt.monte_carlo_runs):
                if self.opt.stop_requested:
                    return []
                dot_progress(n)
                if self.opt.interpoint_relations['simulated - simulated']:
                        distlis = self.__get_same_interpoint_distances(
                            mcli[n]['pli'])
                        mcli[n]['simulated - simulated']['dist'].append(
                            distlis[0])
                        mcli[n]['simulated - simulated']['latdist'].append(
                            distlis[1])
                if self.opt.interpoint_relations['simulated - particle']:
                    distlis = self.__get_interpoint_distances2(mcli[n]['pli'], pli)
                    mcli[n]['simulated - particle']['dist'].append(distlis[0])
                    mcli[n]['simulated - particle']['latdist'].append(distlis[1])
                if self.opt.interpoint_relations['particle - simulated']:
                    distlis = self.__get_interpoint_distances2(pli, mcli[n]['pli'])
                    mcli[n]['particle - simulated']['dist'].append(distlis[0])
                    mcli[n]['particle - simulated']['latdist'].append(distlis[1])
        if self.opt.determine_clusters:
            sys.stdout.write("\n-> Determining simulated clusters:\n   ")
            for n, li in enumerate(mcli):
                if self.opt.stop_requested:
                    return []
                dot_progress(n)
                mcli[n]['clusterli'] = self.__determine_clusters(li['pli'])
                self.__process_clusters(mcli[n]['clusterli'])
        self.mcli = mcli
        sys.stdout.write("\n")
        return mcli

    def __process_clusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.convex_hull = geometry.convex_hull(c)
            hull_centroid = c.convex_hull.centroid()
            c.dist_to_path = hull_centroid.perpend_dist_closed_path(self.path)
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.nearest_cluster = ClusterData()
            if len(clusterli) == 1:
                c.dist_to_nearest_cluster = -1
                return
            c.dist_to_nearest_cluster = sys.maxsize
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateral_dist_to_cluster(c2, self.path)
                    if d < c.dist_to_nearest_cluster:
                        c.dist_to_nearest_cluster = d
                        c.nearest_cluster = c2

    def __determine_clusters(self, pointli):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than opt.within_cluster_dist from at least one
            other point in the cluster
        """
        if self.opt.within_cluster_dist < 0:
            return
        clusterli = []
        for p1 in pointli:
            if self.opt.stop_requested:
                return []
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.to_pixel_units(self.opt.within_cluster_dist,
                                                                       self.pixelwidth):
                    if p2.cluster is not None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        self.__process_clusters(clusterli)
        return clusterli

    def __parse(self):
        """ Parse profile data from input file 
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = file_io.read_file(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        while li:
            s = li.pop(0).replace('\n', '').strip()
            if s.split(' ')[0].upper() == 'IMAGE':
                self.src_img = s.split(' ')[1]
            elif s.split(' ')[0].upper() == 'PROFILE_ID':
                try:
                    self.id = int(s.split(' ')[1])
                except (IndexError, ValueError):
                    profile_warning(self, "Profile id not defined or invalid")
            elif s.split(' ')[0].upper() == 'COMMENT':
                try:
                    self.comment = s.split(' ', 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(' ')[0].upper() == 'PIXELWIDTH':
                try:
                    self.pixelwidth = float(s.split(' ')[1])
                    self.metric_unit = s.split(' ')[2]
                except (IndexError, ValueError):
                    raise ProfileError(self, "PIXELWIDTH is not a valid number")
            elif s.upper() in ("PLASMA_MEMBRANE", "PROFILE_BORDER"):
                self.path = ProfileBorderData(self.__get_coords(li, 'path'))
            elif s.upper() in ("POSTSYNAPTIC_DENSITY", "PSD_OR_ACTIVE_ZONE"):
                self.psd_li.append(Psd(self.__get_coords(li, "PSD"), self))
            elif s.upper() in ("PROFILE_HOLE", "HOLE"):
                self.holeli.append(geometry.SegmentedPath(self.__get_coords(li, 'hole')))
            elif s.upper() in ("POINTS", "PARTICLES"):
                self.pli = PointList(self.__get_coords(li, "particle"), "particle", self)
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(self.__get_coords(li, "random"), "random", self)
            elif s.upper() == 'GRID':
                # Retrieve coordinates to dummy variable as they will not be used
                __ = PointList(self.__get_coords(li, 'grid'), 'grid', self)
                profile_warning(self, "Grid found; however, as grids are no longer supported " 
                                      "it will be discarded")
            elif s[0] != "#":  # unless specifically commented out
                profile_warning(self, "Unrecognized string '" + s +
                                "' in input file")
        # Now, let's see if everything was found
        self.__check_parsed_data()

    def __check_parsed_data(self):
        """See if the profile data was parsed correctly, and print info
        on the parsed data to stdout.
        """
        self.__check_var_default('src_img', "Source image", "N/A")
        self.__check_var_default('id', "Profile ID", "N/A")
        self.__check_var_default('comment', "Comment", "")
        self.__check_var_val('metric_unit', "Metric unit", 'metric_unit')
        self.__check_required_var('pixelwidth', "Pixel width", self.metric_unit)
        # self._check_var_default('profile_type', "Profile type", "N/A")
        self.__check_list_var('path', 'Profile border', 'nodes', 2)
        self.__check_table_var('psd_li', 'Postsynaptic density',
                               'Postsynaptic densities', 1, 2)
        self.__check_list_var('pli', 'Points', '', 0)
        self.__check_table_var('holeli', "Hole", "Holes", 0, 2)
        self.__check_var_exists('randomli', "Random points", 'use_random')

    def __check_required_var(self, var_to_check, var_str, post_str):
        """Confirm that self has a required variable; else, raise
        ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n" % (var_str, self.__dict__[var_to_check], post_str))

    @staticmethod
    def __check_list_len(var, min_len):
        """Return True if var is a list and has at least min_len
        elements, else False.
        """
        return isinstance(var, list) and len(var) >= min_len

    def __check_list_var(self, var_to_check, var_str, post_str, min_len):
        """Confirms that self has a var_to_check that is a list and
        has at least min_len elements; if var_to_check does not exist
        and min_len <= 0, assigns an empty list to var_to_check. Else,
        raise a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file" % var_str)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n" % (var_str, post_str, len(self.__dict__[var_to_check])))

    def __check_table_var(self, var_to_check, var_str_singular,
                          var_str_plural, min_len_1, min_len_2):
        """Confirms that var_to_check exists, is a list and has at
        least min_len_1 elements, and that each of these has at least
        min_len_2 subelements; if var_to_check does not exist and
        min_len_1 <= 0, assigns an empty list to var_to_check. Else,
        raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file" % var_str_plural)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file" % var_str_plural.lower())
        else:
            for element in self.__dict__[var_to_check]:
                if not self.__check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates" % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural, len(self.__dict__[var_to_check])))

    def __check_var_default(self, var_to_check, var_str, default=""):
        """Checks if var_to_check exists; if not, assign the default
        value to var_to_check. Never raises a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            self.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str, self.__dict__[var_to_check]))

    def __check_var_exists(self, var_to_check, var_str, optflag):
        """Checks for consistency between profiles with respect to the
        existence of var_to_check (i.e., var_to_check must be present
        either in all profiles or in none).

        If optflag is not set (i.e., this is the first profile), then
        set optflag to True or False depending on the existence of
        var_to_check. If optflag is already set (for consequent
        profiles), var_to_check must (if optflag is True) or must not
        (if optflag is False) exist. If not so, raise ProfileError.
        """
        if not hasattr(self.opt, optflag):
            if self.__dict__[var_to_check]:
                self.opt.__dict__[optflag] = True
            else:
                self.opt.__dict__[optflag] = False
        if self.opt.__dict__[optflag]:
            if self.__dict__[var_to_check]:
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s expected but not found in input file" % var_str)
        elif self.__dict__[var_to_check]:
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def __check_var_val(self, var_to_check, var_str, optvar):
        """Checks for consistency between profiles with respect to the
        value of var_to_check (i.e., var_to_check must be present and
        have equal value in all profiles).

        If optvar is not set (i.e., this is the first profile), then
        set optflag to the value of var_to_check. If optvar is already
        set (for consequent profiles), the value of var_to_check must
        be equal to that of optvar. If not so, raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(self.opt, optvar):
            self.opt.__dict__[optvar] = self.__dict__[var_to_check]
        elif self.__dict__[var_to_check] == self.opt.__dict__[optvar]:
            pass  # really no point in pointing out that it's ok
            # sys.stdout.write("  %s: %s\n"
            #                  % (var_str, parent.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, self.__dict__[var_to_check],
                                  self.opt.__dict__[optvar]))

    def __check_paths(self):
        """Check if profile border and holes intersect with themselves."""

        def check_path(_path, s):
            for p in range(0, len(_path) - 3):
                for q in range(0, len(_path) - 1):
                    if p not in (q, q + 1) and p + 1 not in (q, q + 1):
                        if geometry.segment_intersection(_path[p],
                                                         _path[p + 1],
                                                         _path[q],
                                                         _path[q + 1]):
                            raise ProfileError(
                                self, "%s invalid (crosses itself)" % s)
            return True

        check_path(self.path, "Profile border")
        for path in self.psd_li:
            check_path(path, "PSD")
        for path in self.holeli:
            check_path(path, "Hole")
        for n, h in enumerate(self.holeli):
            if not h.is_simple_polygon():
                raise ProfileError(self, "Profile hole %d is not a simple polygon" % (n + 1))
            for n2, h2 in enumerate(self.holeli[n + 1:]):
                if h.overlaps_polygon(h2):
                    raise ProfileError(self, "Profile hole %d overlaps with hole %d "
                                       % (n + 1, n + n2 + 2))
        sys.stdout.write("  Paths are ok.\n")

    def get_total_synm(self):
        """Construct a path comprising the plasma membrane of all
        PSDs/active zones and perforations.
        Assume well-behaved PSDs.
        """

        def dist_to_left_endnode(_p):
            path = geometry.SegmentedPath()
            project, seg_project = _p.project_on_path_or_endnode(self.path)
            path.extend([self.path[0], project])
            for n in range(1, seg_project):
                path.insert(len(path)-1, self.path[n])
            return path.length()

        left_mindist = right_mindist = self.path.length()
        left_p = geometry.Point()
        right_p = geometry.Point()
        for p in itertools.chain(*[[psd_[0], psd_[-1]] for psd_ in self.psd_li]):
            left_d = dist_to_left_endnode(p)
            right_d = self.path.length() - left_d
            if left_d <= left_mindist:
                left_mindist = left_d
                left_p = p
            elif right_d <= right_mindist:
                right_mindist = right_d
                right_p = p
        pseudo_psd = Psd(geometry.SegmentedPath([left_p, right_p]), self)
        return pseudo_psd.get_synm()

    def __get_coords(self, strli, coord_type=""):
        """Pop point coordinates from list strli.

        When an element of strli is not a valid point, a warning is
        issued.
        """
        pointli = []
        s = strli.pop(0).replace('\n', '').replace(' ', '').strip()
        while s != 'END':
            try:
                p = geometry.Point(float(s.split(',')[0]), float(s.split(',')[1]))
                if pointli and (p == pointli[-1] or (coord_type == 'particle' and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coord_type, p))
                else:
                    pointli.append(Point(p.x, p.y, ptype=coord_type))
            except ValueError:
                if s[0] != '#':
                    profile_warning(self, "'%s' not valid %s coordinates" % (s, coord_type))
                else:
                    pass
            s = strli.pop(0).replace('\n', '').strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli
# end of class Profile


class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 200
        self.outputs = {'profile summary': True,
                        'particle summary': True,
                        'random summary': True,
                        'session summary': True}
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
        self.use_random = False
        self.determine_clusters = False
        self.within_cluster_dist = 50
        self.run_monte_carlo = False
        self.monte_carlo_runs = 99
        self.monte_carlo_simulation_window = 'profile'
        self.monte_carlo_strict_location = False
        self.determine_interpoint_dists = False
        self.interpoint_dist_mode = 'nearest neighbour'
        self.interpoint_relations = {'particle - particle': True,
                                     'random - particle': True,
                                     'particle - simulated': False,
                                     'simulated - particle': False,
                                     'simulated - simulated': False}
        self.interpoint_shortest_dist = True
        self.interpoint_lateral_dist = False
        self.stop_requested = False

    def reset(self):
        """ Resets all options to default, and removes those that are not
            set in __init__().
        """
        self.__dict__ = {}
        self.__init__()
# end of class OptionData


class ProfileError(Exception):
    def __init__(self, profile, msg):
        self.profile = profile
        self.msg = msg + "."


def profile_warning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True


def profile_message(msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)
