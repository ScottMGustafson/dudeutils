import xml.etree.ElementTree as et

import numpy as np
from scipy import constants

import dudeutils.xmlutils as xmlutils
from dudeutils.atomic import *

tf = {"true": True, "false": False}
c = constants.c / 1000.  # speed of light in km/s


class ObjList(list):
    """
    A container class inheriting from list.
    """

    def __init__(self, objlist, *args, **kwargs):
        super(ObjList, self).__init__()
        self.objlist = objlist
        for key, val in kwargs.items():
            try:
                assert (type(val) in [str, float, int])
            except AssertionError:
                raise AssertionError("invalid type specified: %s:%s" % (
                    str(type(val)), str(val)))
            setattr(self, key, val)

    def __str__(self):
        msg = "%s id:%s\n %s" % (str(type(self)), str(id(self)), str([str(it) for it in self.objlist]))
        return msg

    def __eq__(self, other):
        if not isinstance(other, type(self)) or len(self.objlist) != len(other.objlist):
            return False
        for i in range(len(other.objlist)):
            if other.objlist[i] != self.objlist[i]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, rhs):
        self.objlist = self.objlist + rhs.objlist
        return self

    def __iter__(self):
        for i in range(len(self.objlist)):
            yield self.objlist[i]

    def __getitem__(self, i):
        return self.objlist[i]

    def __delitem__(self, i):
        del self.objlist[i]

    def __contains__(self, item):
        return item in self.objlist

    def __setitem__(self, i, value):
        self.objlist[i] = value

    def __len__(self):
        return len(self.objlist)

    @staticmethod
    def split_types(obj_list):
        types = set([type(it) for it in obj_list])

        ret_lst = []
        for typ in types:
            _lst = [it for it in obj_list if isinstance(it, typ)]
            _new = ObjList.factory(objlist=_lst)
            assert len(_new.objlist) > 0 and _new.objlist != []

            ret_lst.append(_new)

        for item in ret_lst:
            assert len(item.objlist) > 0
            assert item.objlist != [], str(type(item)) + 'len = ' + str(len(item.objlist)) + ' ' + str(item.objlist)
        return ret_lst

    def get_item(self, iden):
        """get element from objlist.
        input:
        ------
        iden: id of elemnet to fetch

        output:
        -------
        data_types.Data subclass instance

        raises:
        -------
        none
        """
        for item in self.objlist:
            if iden == item.id:
                return item

    @staticmethod
    def factory(objlist=None, **kwargs):
        """
        factory method to produce instances of the chidren of ObjList.  

        Input:
        ------
        objlist : a list of data.  elements should be instances of some derived 
                  class of Data

        Output:
        -------
        object instance of desired subclass

        Raises:
        -------
        TypeError when cls isn't recognized

        """
        if objlist == None or objlist == []:
            return None

        if len(set([type(it) for it in objlist])) > 1:
            raise Exception('Mixed types in objlist')

        for cls in ObjList.__subclasses__():
            if cls.registrar_for(ObjList.classname(objlist[0].__class__)):
                return cls(objlist, **kwargs)
        raise TypeError("invalid type: " + ObjList.classname(objlist[0].__class__))

    @staticmethod
    def classname(cls):
        return cls.__name__.split('.')[-1]

    def xml_rep(self, parent):
        """return the list of all relevant nodes in xml"""
        if not isinstance(parent, str):
            parent = ObjList.classname(
                parent.__class__)  # if a class is given instead of class name, just get the class name
            parent = parent.split('.')[-1]
        current = et.SubElement(parent, ObjList.classname(self.__class__), {"id": self.id})
        current.extend([item.node for item in self.objlist])
        return current

    @staticmethod
    def get_class(node):
        """
        get class of specified node.

        Input:
        ------
        node : xml node of data

        Output:
        -------
        the class that matches the specified node tag

        Raises:
        -------
        Exception if tag not recognized
        """
        for item in ObjList.__subclasses__():
            if node.tag == ObjList.classname(item):
                return item
        raise Exception("class %s not recognized" % (node.tag))

    @staticmethod
    def _read_node(node, verbose=False):
        """
        given a node that is an AbsorberList or similar type parse individual 
        datum

        input:
        ------
        node: node of one of AbsorberLists, etc...

        output:
        -------
        ObjList (or subclass) instance
        """
        if verbose:
            print(len(list(node)))
        objlist = [Data.factory(node=item) for item in list(node)]  # each item should be xml node
        return ObjList.factory(objlist, id=node.get("id"))

    @staticmethod
    def sublass_str():
        """return subclass names as list of strings"""
        return [item.__name__ for item in ObjList.__subclasses__()]

    @staticmethod
    def get_from_xml(id, parent):
        """
        returns ObjList object instance from a given parent.  

        Input:
        ------
        parent : xml absorberLists node
        id : unique identifier of the data.

        Output:
        -------
        ObjList instance 

        """
        for item in list(parent):
            if item.tag in [it + "s" for it in ObjList.subclass_str()]:
                if item.id == id:
                    return ObjList._read_node(item)
        return None

    def append_datum(self, tag, **kwargs):
        node = et.Element(tag, **kwargs)
        self.objlist.append(Data.factory(node=node))

    def refresh_list(xmlfile):
        tree = et.parse(xmlfile)
        root = tree.getroot()
        for tag in "AbsorberLists ContinuumPointLists RegionLists SingleViewLists VelocityViewLists".split():
            parent = root.find(tag)
            for objlist in parent:
                inst = ObjList._read_node(objlist)

    @staticmethod
    def list_from_xml(parent, verbose=False):
        """
        returns a list of ObjList objects from a given parent. 

        Input:
        ------
        parent : xml absorberLists node (list of absorberList instances)

        Output:
        -------
        list of class instances derived from ObjList 

        """
        if verbose:
            print(len(parent), ' elements')
        return [ObjList._read_node(item, False) for item in parent]


class AbsorberList(ObjList):
    @classmethod
    def registrar_for(cls, tag):
        return tag == "Absorber"

    def patch_list(self, lst):
        """suppose you have some list with some absorbers to replace..  then use this function, my friend."""
        while len(lst) > 0:
            new_val = lst.pop(0)
            if not new_val.id:
                continue
            for i in range(len(self.objlist)):
                if self.objlist[i].id == new_val.id:
                    self.objlist[i] = new_val


class ContinuumPointList(ObjList):
    @classmethod
    def registrar_for(cls, tag):
        return tag == "ContinuumPoint"

    def patch_list(self, lst):
        """suppose you have some list with some absorbers to replace..  then use this function, my friend."""
        while len(lst) > 0:
            new_val = lst.pop(0)
            if not new_val.id:
                continue
            for i in range(len(self.objlist)):
                if self.objlist[i].id == new_val.id:
                    self.objlist[i] = new_val


class RegionList(ObjList):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.consolidate_regions()

    def in_regions(self, val):
        for reg in self.objlist:
            if reg.start <= val <= reg.end:
                return True
        return False

    @classmethod
    def registrar_for(cls, tag):
        return tag == "Region"

    @staticmethod
    def consolidate_list(ranges):
        """
        consolidate list of ranges to combine overlapping ranges as one.

        Input:
        ------
        ranges:  list of ranges  (list of length=2 lists of floats)

        Output:
        -------
        result: list of consolidated ranges  (list of length=2 lists of floats)

        Raises:
        -------
        None

        """
        try:
            assert isinstance(ranges, list)
        except:
            raise TypeError('argument should be list-like.  Instead got %s' % str(type(ranges)))
        result = []
        current_start = -1
        current_stop = -1

        for start, stop in sorted(ranges):
            if start > current_stop:
                # this segment starts after the last segment stops
                # just add a new segment
                result.append((start, stop))
                current_start, current_stop = start, stop
            else:
                # segments overlap, replace
                result[-1] = (current_start, stop)
                # current_start already guaranteed to be lower
                current_stop = max(current_stop, stop)
        return result

    def consolidate_regions(self):
        """
        adapted from:
        http://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap
        """
        self.objlist = sorted(self.objlist, key=lambda x: x.start)
        ranges = [(item.start, item.end) for item in self.objlist]

        result = RegionList.consolidate_list(ranges)

        for i in range(len(self.objlist)):
            if i < len(ranges):
                self.objlist[i].start = ranges[i][0]
                self.objlist[i].end = ranges[i][-1]

        while len(self.objlist) > len(ranges):
            del (self.objlist[-1])

    def get_indices(self, waves):
        indices = []
        for region in self.objlist:
            ind1 = np.where(waves < region.end)[0]
            ind2 = np.where(waves > region.start)[0]
            indices += list(set(ind1).intersection(ind2))
        return sorted(list(set(indices)))


class SingleViewList(ObjList):
    @classmethod
    def registrar_for(cls, tag):
        return tag == "SingleView"


class VelocityViewList(ObjList):
    @classmethod
    def registrar_for(cls, tag):
        return tag == "VelocityView"


class Data(object):
    node_attrib = []

    def __init__(self, *args, **kwargs):
        tag = kwargs.get('tag', False)
        if not tag:
            self.tag = self.__class__.__name__
        self.node = kwargs.pop('node', None)
        if self.node:
            self.parse_node()
        else:
            for key, val in kwargs.items():
                try:
                    setattr(self, key, float(val))
                except (ValueError, TypeError) as e:
                    if val in ['true', 'false']:
                        val = True if val == 'true' else False
                    setattr(self, key, val)

                    # self.get_mock_node()

    def get_mock_node(self):
        attrs = [attr for attr in dir(self) if not attr.startswith('__') and attr not in ['tag', 'node']]
        attrib = {str(k): str(getattr(self, k)) for k in attrs}
        self.node = type('mock_node', (), {"tag": 'Absorber', 'attrib': attrib})

    def validate_init(self):
        pass

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        return self.node.attrib == other.node.attrib

    def __ne__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def factory(**kwargs):
        tag = kwargs.get("tag", None)
        if not tag:
            try:
                tag = kwargs.get("node").tag
                assert tag is not None, 'no tag found.  Is this a valid xml node? %s' % str(kwargs.get("node", None))
            except:
                raise Exception("need to specify either tag and id or node")

        for cls in Data.__subclasses__():
            if cls.registrar_for(tag):
                if "node" in kwargs.keys():
                    inst = cls.from_node(**kwargs)
                elif "xmlfile" in kwargs.keys():
                    inst = cls.from_file(**kwargs)
                else:
                    return cls(**kwargs)

                # inst.keys = list(inst.node.attrib.keys())
                # inst.parse_node()
                try:  # if an additional constructor is specified, run it
                    inst.alt_init(**kwargs)
                except:
                    pass
                return inst

        raise ValueError("\n%s not a valid data type.  \n"
                         "Valid types are \n  %s" % (str(tag), str(Data.__subclasses__())))

    @classmethod
    def from_file(cls, **kwargs):
        """constructor from file"""
        tag = kwargs.pop("tag")
        id = kwargs.pop("id")
        xmlfile = kwargs.pop("xmlfile")

        root = et.parse(xmlfile).getroot()
        node = xmlutils.get_node(root.find('CompositeSpectrum'), tag, id)

        return cls(tag, id=id, xmlfile=xmlfile, node=node, **kwargs)

    @classmethod
    def from_node(cls, **kwargs):
        """constructor from node"""
        node = kwargs.get("node")
        inst = cls(**kwargs)
        inst.parse_node()
        inst.validate_init()

        return inst

    @staticmethod
    def get_node_attrib(cls, kwargs):
        _kwargs = kwargs
        for key in kwargs.keys():
            if not key in cls.node_attrib:
                del (_kwargs[key])
        return _kwargs

    @staticmethod
    def read(xml_data, tag=None, ids=None):

        """
        read the data.
        xml_data: input xml_data, either file-like object, filename or string of xml data
        tag (default='Absorber'): the specified data type
        ids: which ids to read.  if None, then read all available.
        """

        if xml_data.endswith('.xml'):
            the_tree = et.parse(xml_data)
        else:
            assert isinstance(xml_data, str)
            the_tree = et.ElementTree(et.fromstring(xml_data))

        duderoot = the_tree.getroot()  # should be SpecTool

        dudespec = duderoot.find("CompositeSpectrum")
        if dudespec is None:
            raise Exception("error reading fit file.  check %s to verify." % (xml_data))
        # spectrum = dudespec.find("Spectrum")
        lst = []
        specdata = list(dudespec) + list(duderoot)  # + list(the_tree.findall('Region'))

        if ids is None:
            for item in specdata:
                if item in lst:
                    continue
                if tag:
                    if item.tag == tag:
                        lst.append(Data.factory(node=item))
                else:
                    if item.tag in "Region Absorber ContinuumPoint".split(" "):
                        lst.append(Data.factory(node=item))

        else:
            try:
                if len(ids) == 0:
                    raise Exception('need at least one id specified')
            except:
                raise Exception("type of ids needs to be list of str")
            for item in specdata:
                if item in lst:
                    continue
                if item.tag == tag and item.get('id') in ids:
                    lst.append(Data.factory(node=item))
        return lst

    def locked(self, param):
        return getattr(self, param + "Locked")

    def parse_node(self, node=None):
        """read from node, set attribs to self"""
        if node == None:
            node = self.node

        data = node.attrib
        for key, val in data.items():
            assert isinstance(key, str), 'got: %s' % str(key)
            if "Locked" in key:
                if isinstance(val, str):
                    setattr(self, key, tf[val.lower()])
            else:
                try:
                    setattr(self, str(key), float(val))
                except:
                    setattr(self, str(key), str(val))

    def set_node(self, **kwargs):
        """set values from self to node"""
        for key, val in dict(kwargs).items():
            self.node.set(key, str(val))

    def set_data(self, **kwargs):
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()):
            if "Locked" in key:
                if not type(val) is bool:
                    assert isinstance(val, str)
                    assert val in tf.keys()
                self.node.set(key, str(val).lower())
            setattr(self, key, val)
            self.node.set(key, str(val))


class Absorber(Data):
    node_attrib = ["id", "ionName",
                   "N", "NLocked", "NError",
                   "b", "bLocked", "bError",
                   "z", "zLocked", "zError"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.defaults()

    def __eq__(self, other):
        if not isinstance(other, Absorber):
            return False
        for attr in "id ionName N b z".split(" "):
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def validate_init(self):
        for attr in Absorber.node_attrib:
            if not hasattr(self, attr):
                if attr in "NLocked bLocked zLocked".split(' '):
                    setattr(self, attr, True)
                elif attr == 'id':
                    setattr(self, attr, '')
                elif attr in "NError bError zError".split(' '):
                    setattr(self, attr, 0.0)
                else:
                    raise Exception("Absorber instance requires setting of ionName, N,b,z")

    def defaults(self):
        for attr in Absorber.node_attrib:
            if not hasattr(self, attr):
                if attr in "NLocked bLocked zLocked":
                    setattr(self, attr, True)
                elif attr == 'id':
                    setattr(self, attr, '')
                elif attr in "NError bError zError":
                    setattr(self, attr, 0.0)

    @classmethod
    def registrar_for(cls, tag):
        return tag == "Absorber"

    def __str__(self):
        return "%-5s N=%8.5lf b=%8.5lf z=%10.8lf, NLocked=%s, bLocked=%s, zLocked=%s" % (
            self.ionName, self.N, self.b, self.z, str(self.NLocked), str(self.bLocked), str(self.zLocked))

    def getShift(self, z):
        return (float(self.z) - z) * c / (1. + z)

    def get_wave(self, n=0):
        """get observed wave.  n=transition level such that 
            lya (n=2-->n=1)=0
            lyb (n=3-->n=1)=1
            and etc..
        """
        return (1. + self.z) * atomic_data[self.ionName][n].wave

    def get_obs_wave(self, n=0):
        return self.get_wave(n)

    def get_rest_wave(self, n=0):
        """get observed wave.  n=transition level such that 
            lya (n=2-->n=1)=0
            lyb (n=3-->n=1)=1
            and etc..
        """
        return atomic_data[self.ionName][n].wave

    def get_lines(self):
        """
        return list of SpectralLine instances containing relevant atomic data
        """
        lst = []
        if not self.ionName in atomic_data.keys():
            print(self.ionName + ": not in keys")
            print("available keys are: \n" + str(atomic_data.keys()))
            raise KeyError()
        for item in atomic_data[self.ionName]:
            kwargs = {}
            for key in ['f', 'gamma']:
                kwargs[key] = getattr(item, key)
            for key in Absorber.node_attrib:
                try:
                    kwargs[key] = getattr(self, key)
                except:
                    kwargs[key] = ''
            kwargs['wave'] = float(getattr(item, 'wave'))
            kwargs['obs_wave'] = (1. + self.z) * float(getattr(item, 'wave'))
            lst.append(SpectralLine(**kwargs))
        return lst

    # @staticmethod
    # def get_lyman_limit_tau(wave, z, N):
    #     """get tau for lyman limit systems"""
    #     lambda_l = 911.8
    #     if wave > lambda_l - 0.1:
    #         wave = lambda_l - 0.1;
    #
    #     acotz = np.atan(1. / z)
    #     t2 = np.exp(-4. * z * acotz) / (1. - np.exp(-2. * np.pi * z))
    #     gaunt = 8. * np.pi * np.sqrt(3.) * (wave / lambda_l) * t2;
    #     return np.pow(10., N) * ((wave / lambda_l) ** 3.) * 7.91e-18 * gaunt;

    def locked(self, param):
        param_lock = {'N': 'NLocked', 'b': 'bLocked', 'z': 'zLocked'}
        try:
            ans = getattr(self, param_lock[param])
            if str(ans) not in ['true', 'True', 'TRUE']:
                return False
            else:
                return True
        except KeyError:
            raise Exception("no param named %s" % {param})

    def in_region(self, regions):
        """
        detects whether of not CENTER of line is in a optimization region.
        only needs to be in one region to pass
        """
        for region in regions:
            for line in self.get_lines():
                if line.get_obs(z=line.z) in region and line.get_equiv_width(line.N, line.z, pixels=True) > 0.29:
                    return True
        return False


class ContinuumPoint(Data):
    node_attrib = ["id", "x", "xLocked", "xError", "y", "yLocked", "yError"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.validate_init()
        self.defaults()

    def __eq__(self, other):
        if not isinstance(other, ContinuumPoint):
            return False
        for attr in ContinuumPoint.node_attrib:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def validate_init(self):
        for attr in ContinuumPoint.node_attrib:
            if not hasattr(self, attr):
                if 'Locked' in attr:
                    setattr(self, attr, True)
                elif attr == 'id':
                    setattr(self, attr, 'null')
                elif 'Error' in attr:
                    setattr(self, attr, 0.0)

    def defaults(self):
        for attr in ContinuumPoint.node_attrib:
            if not hasattr(self, attr):
                if attr in "xLocked yLocked":
                    setattr(self, attr, True)
                elif attr == 'id':
                    setattr(self, attr, 'null')
                elif attr in "xError yError":
                    setattr(self, attr, 0.0)

    def __str__(self):
        return "id=%10s x=%7.2lf y=%10.8E xLocked=%s yLocked=%s" % (
            self.id, self.x, self.y, str(self.xLocked), str(self.yLocked))

    @classmethod
    def registrar_for(cls, tag):
        return tag == "ContinuumPoint"

    def in_region(self, regions):
        for region in regions:
            if self.x in region:
                return True
        return False


class Region(Data):
    node_attrib = ["start", "end"]

    def __init__(self, *args, **kwargs):
        self.start = kwargs.get('start')
        self.end = kwargs.get('end')
        super().__init__(*args, **kwargs)

    def validate_init(self):
        for attr in Region.node_attrib:
            if not hasattr(self, attr):
                raise Exception("Region instance requires setting of start, end")

    def __eq__(self, other):
        if not isinstance(other, Region):
            return False
        return self.start == other.start and self.end == other.end

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, val):
        return self.start <= float(val) <= self.end

    def __gt__(self, other):
        return self.end > other.end

    def __lt__(self, other):
        return self.start < other.start

    def __ge__(self, other):
        return self.end >= other.end

    def __le__(self, other):
        return self.start <= other.start

    def __str__(self):
        return "start=%lf, end=%lf" % (self.start, self.end)

    @classmethod
    def registrar_for(cls, tag):
        return tag == "Region"


class SingleView(Data):
    node_attrib = ["id", "centWave", "waveRange", "minFlux", "maxFlux"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def registrar_for(cls, tag):
        return tag == "SingleView"


class VelocityView(Data):
    node_attrib = ["id", "labels",
                   "minWave", "maxWave", "minFlux", "maxFlux",
                   "restWaves", "redshift"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def registrar_for(cls, tag):
        return tag == "VelocityView"


class Param(object):
    """
    class for the manipulation of individual absorption parameters, N,b and z"""

    def __init__(self, param_name, value, locked, error=0., parent=None, index=None, bounds=None):
        self.name = param_name
        self.absorber = parent
        self.locked = locked
        self.error = float(error)
        self.value = float(value)
        self.index = int(index)
        if bounds:
            self.guess = Param.random_initial_cond(bounds)
        else:
            self.guess = float(value)

    @staticmethod
    def random_initial_cond(bounds):
        """
        provide randomized initial conditions given some set of bounds

        Input:
        ------
        bounds: length=2 list of floats

        Output:
        -------
        (float)

        Raises:
        -------
        None
        """
        return (bounds[1] - bounds[0]) * np.random.random_sample() + bounds[0]

    def __str__(self):
        return "%s=%lf, %sError=%lf, %sLocked=%s" % (
            self.name, self.value,
            self.name, self.error,
            self.name, str(self.locked).lower())


if __name__ == '__main__':
    print('instantiating an ab')
    ab = Absorber(id='id', ionName='C IV',
                  N=13.45, NLocked=True, NError=0,
                  b=13.45, bLocked=True, bError=0,
                  z=3.45)
    print('cool.  instantiating another with a node.')
    attrib_dct = {'ionName': 'C IV', 'N': float(11.2), 'b': float(2.02), 'z': float(1.92), 'NError': 0, 'bError': 0,
                  'zError': 0.0, 'NLocked': False, 'bLocked': False, 'zLocked': False, 'id': ""}
    test_node = type('testnode', (), {"tag": 'Absorber', 'attrib': attrib_dct})
    ab = Absorber(node=test_node)
