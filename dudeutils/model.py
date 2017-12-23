import copy
import os.path
import pickle
import warnings
import xml.etree.ElementTree as et

import numpy as np
from numpy.random import random_sample, randn

import dudeutils.data_types as data_types
from dudeutils.constraints import Constraint
from dudeutils.wavelength import get_waves

c = 299792.458  # speed of light in km/s


def tf(val):
    return "true" if val in [True, "true"] else "false"


class Model(object):
    # this dict maps ObjList subclasses to their associated data type
    model_classes = {"AbsorberList": "Absorber",
                     "ContinuumPointList": "ContinuumPoint",
                     "RegionList": "Region",
                     "SingleViewList": "SingleView",
                     "VelocityViewList": "VelocityView"}

    def __init__(self, **kwargs):
        self.pixels = 0
        self.chi2 = 0.
        self.params = 0
        self.flux = None  # fits (or text) file with flux
        self.error = None  # fits (or text) file with flux.  if test, then this is the same as flux and xml format will change a bit
        self.xmlfile = kwargs.pop("xmlfile", None)
        self.get_all = kwargs.pop("get_all", True)
        self.abs_ids = kwargs.pop("abs_ids", None)

        if self.xmlfile:
            self.read()
        else:  # lets you instantiate without specifying xml
            self.xmlfile = 'test_file.xml'
            if not "AbsorberList" in kwargs.keys():
                raise Exception("AbsorberList should be specified on instantiation")

        for key, val in kwargs.items():
            setattr(self, key, val)

        self.chi2 = float(self.chi2)

        self.pixels = int(float(self.pixels))
        self.params = int(float(self.params))

        self.locked = {}
        self._dof = float(self.pixels) - float(self.params)
        for k in Model.model_classes.keys():
            if not hasattr(self, k):
                setattr(self, k, [])

    def __eq__(self, other):
        attrs = list(['xmlfile', 'chi2', 'pixels', 'params']) + list(Model.model_classes.keys())
        for item in attrs:
            if getattr(self, item) != getattr(other, item):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        """
        string output for a model will be like:
    
        id=HI N=17.12345 b=12.345678 z=1.234567890
        id=SiI N=11.12345 b=12.345678 z=1.234567890
        id=OI N=11.12345 b=12.345678 z=1.234567890
        locked=OI:bLocked HI:zLocked
        chi2=123.4 pixels=154 params=23

        """
        string = "model : src, flux = %s, %s\n" % (str(self.xmlfile), str(self.flux))
        string += "\nchi2=%lf pixels=%lf params=%lf\n\n" % (
            float(self.chi2), float(self.pixels), float(self.params))
        string += "-----------AbsorberList------------\n"
        for item in self.absorber_list:
            string += str(item) + "\n"

        string += "-----------ContinuumPointList------\n"
        for item in self.cont_point_list:
            string += str(item) + "\n"

        string += "-----------RegionList--------------\n"
        for item in self.region_list:
            string += str(item) + "\n"

        return string

    def get_indices(self):
        waves = self.get_waves()
        self.indices = self.region_list.get_indices(waves)
        return self.indices

    def get_waves(self):
        if self.flux.endswith('.fits'):
            waves = get_waves(self.flux)
        else:
            col1, col2, *cols = np.loadtxt(self.flux, unpack=True)
            waves = col2 if np.all(col1) == 0 else col1
        return waves

    def _update_pixels(self):

        waves = self.get_waves()
        self.indices = self.region_list.get_indices(waves)
        self.pixels = len(self.indices)

    def _update_params(self):
        self.params = 0
        for ab in self.absorber_list:
            if ab.NLocked and ab.bLocked and ab.zLocked:
                continue
            for line in ab.get_lines():
                if self.region_list.in_regions(line.get_obs(ab.z)):
                    # this seems to be the limit that dude uses
                    if line.get_equiv_width(ab.N, ab.z, pixels=True) > 0.29:
                        self.params += [ab.NLocked, ab.bLocked, ab.zLocked].count(False)
                        break  # to prevent multi-counting of same absorber

    def update_dof(self):
        self._update_params()
        self._update_pixels()

    @property
    def reduced_chi2(self):
        self._reduced_chi2 = float(self.chi2) / (self.dof)
        return self._reduced_chi2

    @reduced_chi2.setter
    def reduced_chi2(self, value):
        self._reduced_chi2 = float(value)

    @property
    def dof(self):
        self._dof = float(self.pixels - self.params)
        return self._dof

    @dof.setter
    def dof(self, value):
        warnings.warn("setting dof without adjusting pixels or params")
        self._dof = value

    @property
    def dh(self):
        d = self.get_datum('D', 'Absorber', param='N')
        h = self.get_datum('H', 'Absorber', param='N')
        self._dh = float(d) - float(h)
        return self._dh

    @dh.setter
    def dh(self, value):
        """shouldnt need this"""
        try:
            self._dh = value
        except AttributeError:
            raise Exception("no availble D/H")

    def toggle_locks(self, dct, locked, tag='Absorber'):
        """ 
        toggle locks in model

        inputs:
        -------
        dct:  dict of idens and params:  {key_name:param_to_toggle, ... }
        locked: bool True to lock
        tag: name of object class

        output:
        -------
        None

        raises:
        -------
        None
        """
        for iden, params in dct.items():
            for param in params:
                if not 'Locked' in param:
                    param += 'Locked'
                self.set_val(iden, tag, **{param: bool(locked)})

    def toggle_cont_lock(self, lst, param='y', locked=True):
        lst = dict(zip([item for item in lst], [param for item in lst]))
        self.toggle_locks(lst, locked, tag='ContinuumPoint')

    def lock_all_cont(self, tf=True):
        for item in self.cont_point_list:
            for param in 'xLocked yLocked'.split():
                self.set_val(item, "ContinuumPoint", **{param: tf})

    def append_datum(self, tag, **kwargs):
        data_types.ObjList.append_datum(tag, **kwargs)

    def get_vel(self, id1, id2):
        """
        get velocity shift between two absorbers

        input:
        ------
        id1, id2:  (string): ids of absorbers to use.  z2 will be our reference.
                    if vel>0, then id1 is red of id2
        output:
        -------
        velocity in km/s

        """
        z1 = float(self.get_datum(id1, 'Absorber', 'z'))
        z2 = float(self.get_datum(id2, 'Absorber', 'z'))
        return (z1 - z2) * c / (1. + z2)

    def get_shift(self, id1, id2):
        """
        an alias for get_vel(id1, id2)

        """
        return self.get_vel(id1, id2)

    def copy(self):
        mod = Model(pixels=self.pixels, chi2=self.chi2, params=self.params,
                    flux=copy.copy(self.flux), error=copy.copy(self.error),
                    abs_ids=copy.copy(self.abs_ids),
                    AbsorberList=copy.deepcopy(self.absorber_list),
                    ContinuumPointList=copy.deepcopy(self.cont_point_list),
                    RegionList=copy.deepcopy(self.region_list))

        mod.xmlfile = copy.copy(self.xmlfile)  # this needs to be done outside of instantiation.

        return mod

    def set_absorbers(self, ab_lst):
        for ab in ab_lst:
            for i, ab_ in enumerate(self.absorber_list):
                if ab_.id == ab.id:
                    self.absorber_list[i] = ab
                    break

    def set_cont_points(self, cont_lst):
        for cnt in cont_lst:
            for i, cnt_ in enumerate(self.cont_point_list):
                if cnt_.id == cnt.id:
                    self.cont_point_list[i] = cnt
                    break

    def build_xml(self, raw_data='', spname='', sptype=''):
        """build a dude-style xml for this model"""
        if raw_data == '':
            raw_data = self.flux
            if not raw_data:
                raise Exception("must specify location of fits source data")
        if spname == '':
            spname = raw_data
        sptype = 'fits' if spname.endswith('.fits') else 'ascii1'
        duderoot = et.Element("SpecTool", {"version": "1.0"})
        dudespec = et.SubElement(duderoot, "CompositeSpectrum", {"id": raw_data})
        spec = et.SubElement(dudespec, "Spectrum", {"spec": spname, "spectype": sptype})

        # add spectrum stuff
        for attr in ["ContinuumPointList", "AbsorberList"]:
            try:
                lst = self.get_lst(attr)
            except:
                if getattr(self, attr) is None:
                    raise Exception("no valid continuum points or absorbers for %s" % (self.xmlfile))
                else:
                    raise
            if lst is None:
                raise Exception(attr + " not found for " + self.xmlfile)
            dudespec.extend([item.node for item in lst])

        # the view stuff
        for attr in ["SingleViewList", "VelocityViewList", "RegionList"]:
            try:
                duderoot.extend([item.node for item in self.get_lst(attr)])
            except:
                pass
        return duderoot

    def get_datum(self, iden, tag, param=None):
        """get an individual datum from model's data list."""
        if tag in Model.model_classes.values():
            tag = inv_dict(tag)

        tst = self.get_lst(tag)
        assert tst

        for item in self.get_lst(tag):
            if iden == item.id:
                if param:
                    if "Locked" in param:
                        return bool(getattr(item, param))
                    try:
                        return float(getattr(item, param))
                    except:
                        return getattr(item, param)
                else:
                    return item
        raise Exception("item not found: %s" % (iden))

    def check_vals(self):
        unphysical = {"b": [0.1, 120.], "N": [8.00, 25.00]}

        abslist = self.absorber_list
        for item in abslist:
            assert (isinstance(item, data_types.Absorber))
            for key, val in unphysical.items():
                try:
                    assert (val[0] <= float(getattr(item, key)) <= val[1])
                except:
                    raise Exception(
                        "%s %s: %s has unphysical value of %lf" % (item.ionName, item.id, key, (getattr(item, key))))

    def monte_carlo_set(self, iden, tag, val_range, param, gaussian=False):
        """set a param for Data `id` to a random value in val_range
        params:
        -------
        iden:  id or object instance of the item to be changed
        tag:   absorber, continuum point, etc..  see data_types.py
        val_range:  list or tuple of [min, max].  if gaussian=true, this will be 95% limits
        gaussian:  return random value on gaussian distribution?
        param: which param to alter?

        """
        a, b = min(val_range), max(val_range)
        if gaussian:
            twosigma = (b - a) / 2.
            sigma = twosigma / 2.
            new = sigma * randn() + (a + b) / 2.
        else:
            new = (b - a) * random_sample() + a
        self.set_val(iden, tag, **{param: new})

    def read(self, xml_data=None):
        """read from xml fit file, apply attribs to self

        input:
        ------
        xml_data (default None)  xml_data to read.  if None, reads self.xmlfile

        output:
        -------
        None

        raises:
        -------
        Exception:  when xml parsing fails


        """

        if not xml_data:
            xml_data = self.xmlfile

        # read through the data:



        lst = data_types.Data.read(xml_data)
        objlst_lst = data_types.ObjList.split_types(lst)
        # todo now does not filter out by id
        for it in objlst_lst:
            cls = it[0].__class__.__name__
            lst_cls = cls + 'List'
            setattr(self, lst_cls, it)

        # read thorugh xml file for source data files
        try:
            if not type(xml_data) is str:
                xml_data.seek(0)
            if xml_data.endswith('.xml'):
                duderoot = et.parse(xml_data).getroot()  # should be SpecTool
            else:
                duderoot = et.ElementTree(et.fromstring(xml_data)).getroot()
        except:
            if type(xml_data) is str:
                raise Exception(xml_data + " failed to parse.")
            else:
                raise
        composite_spec = duderoot.find("CompositeSpectrum")
        path = os.path.split(composite_spec.get("id"))[0]

        # find the source data
        spectrum = composite_spec.find("Spectrum")

        try:
            self.flux = os.path.join(path, spectrum.get('spec'))
        except:
            self.flux = composite_spec.get('id')  # id tag will already have full path

        try:
            self.error = os.path.join(path, spectrum.get('error'))
        except:
            self.error = self.flux.replace('.fits', '_e.fits')

        try:
            self.chi2 = float(spectrum.get("chi2"))
            self.params = float(spectrum.get("params"))
            self.pixels = float(spectrum.get("pixels"))
        except:
            self.chi2 = 0.
            self.params = 0.
            self.pixels = 0.

    def set_val(self, iden, tag="Absorber", **kwargs):
        """
        set value of item in model

        input:
        ------
        iden (unspecified type):  the item or the id of the item to set
        tag (string, optional):  the type of iden

        output:
        -------
        None

        raises:
        -------
        TypeError when type of iden or tag is not recognized

        """

        if type(iden) == str:
            if tag in Model.model_classes.values():
                tag = inv_dict(tag)
            else:
                raise TypeError("model.Model.set_val(): unrecognized type: %s" % (tag))

            # get the model's absorber list
            ab_lst = getattr(self, tag)

            for item in ab_lst:  # for item in abslist
                if item.id == iden:  # if the id of the absorber matches
                    item.set_data(**kwargs)

        else:
            tag = iden.__class__.__name__.split('.')[-1]
            if tag in Model.model_classes.values():
                iden.set_data(**kwargs)
            elif tag in list(Model.model_classes.keys()):
                for key, val in dict(kwargs).items:
                    setattr(iden, key, val)
            else:
                raise TypeError("model.Model.set_val(): unrecognized type: %s" % (tag))

    def write(self, filename=None):
        if filename is None:
            filename = self.xmlfile

        with open(filename, "w") as f:
            f.write("<?xml version=\"1.0\"?>\n")
            f.write("<SpecTool version=\"1.0\">\n")
            if self.flux.endswith('.fits'):
                f.write(
                    "<CompositeSpectrum id=\"%s\"><Spectrum spec=\"%s\" error=\"%s\" chi2=\"%lf\" pixels=\"%lf\" params=\"%d\"/>\n" % (
                        self.flux, os.path.split(self.flux)[-1], os.path.split(self.error)[-1], self.chi2, self.pixels,
                        self.params
                    ))
            else:
                f.write(
                    "<CompositeSpectrum id=\"%s\"><Spectrum spec=\"%s\" spectype=\"ascii1\" chi2=\"%lf\" pixels=\"%lf\" params=\"%d\"/>\n" % (
                        self.flux, os.path.split(self.flux)[-1], self.chi2, self.pixels, self.params
                    ))
            for item in self.cont_point_list:
                if not item.id:
                    item.id = "null"
                f.write(
                    "<ContinuumPoint x=\"%lf\" y=\"%E\" xError=\"0.0\" yError=\"0.0\" xLocked=\"%s\" yLocked=\"%s\" id=\"%s\"/>\n" % (
                        item.x, item.y, tf(item.xLocked), tf(item.yLocked), item.id
                    ))

            for item in self.absorber_list:
                f.write(
                    "<Absorber ionName=\"%s\" N=\"%lf\" b=\"%lf\" z=\"%lf\" NError=\"0.0\" bError=\"0.0\" zError=\"0.0\" NLocked=\"%s\" bLocked=\"%s\" zLocked=\"%s\" id=\"%s\"/>\n" % (
                        item.ionName, item.N, item.b, item.z, tf(item.NLocked), tf(item.bLocked), tf(item.zLocked),
                        item.id
                    ))
            f.write("</CompositeSpectrum>\n")

            for item in self.region_list:
                f.write("<Region start=\"%lf\" end=\"%lf\"/>\n" % (item.start, item.end))
            f.write("</SpecTool>")

    @property
    def region_list(self):
        return self.RegionList

    @region_list.setter
    def region_list(self, val):
        try:
            self.RegionList = data_types.ObjList.factory(val, tag='Region')
        except AttributeError:
            raise Exception('can\'t set empty object list')

    @property
    def absorber_list(self):
        return self.AbsorberList

    @absorber_list.setter
    def absorber_list(self, val):
        try:
            self.AbsorberList = data_types.ObjList.factory(val, tag='Absorber')
        except AttributeError:
            raise Exception('can\'t set empty object list')

    @property
    def cont_point_list(self):
        return self.ContinuumPointList

    @cont_point_list.setter
    def cont_point_list(self, val):
        try:
            self.ContinuumPointList = data_types.ObjList.factory(val, tag='ContinuumPoint')
        except AttributeError:
            raise Exception('can\'t set empty object list')

    def get_spectral_line(self, iden, transition):
        if type(iden) is int:  # gave an index instead of id
            ab = self.absorber_list[iden]
        else:
            ab = self.absorber_list.get_item(iden)
        try:
            return ab.get_lines()[transition]
        except AttributeError:
            raise Exception("absorber %s not found" % (str(iden)))
        except KeyError:
            raise KeyError("absorber %s has no transition %d" % (str(iden), transition))

    @staticmethod
    def _get_item(_id, lst):
        for i, it in enumerate(lst):
            if it.id == _id:
                return i
        raise Exception(_id, ' not in list ', [str(it) for it in lst])

    def get_ab_index(self, _id):
        return Model._get_item(_id, self.absorber_list)

    def get_cont_index(self, _id):
        return Model._get_item(_id, self.cont_point_list)


class ModelDB(list):
    def __init__(self, name=None, models=[], constraints=None, **kwargs):
        """
        Model Database

        Inputs:
        -------
        models:  list of Model instances
        constraints: see ModelDB.constrain
        name: name of the xml models file.  (not the fit file)
        """
        super(ModelDB, self).__init__()

        for key, val in dict(kwargs).items():
            setattr(self, key, val)

        self.name = name
        # self.dbxml=xmlutils.Model_xml(filename=name)
        if len(models) > 0:  # instantiate new db from models
            self.models = models
            for key in list(Model.model_classes.keys()):
                try:
                    setattr(self, key, [getattr(item, key) for item in self.models])
                except:
                    setattr(self, key, [])

        elif name:
            if name.endswith('.xml'):
                self.models = ModelDB.read(str(name), returndb=False)


            else:
                self.models = ModelDB.load_models(name).models
        else:
            self.models = []

        if constraints:
            self.constrain(constraints)

    def __iter__(self):
        for i in range(len(self.models)):
            yield self.models[i]

    def __len__(self):
        return len(self.models)

    def __getitem__(self, i):
        return self.models[i]

    def __eq__(self, other):
        """doensn't matter whether or not other is list or ObjList"""
        for item in other:
            if not item in self.models:
                return False
        return len(other) == len(self.models)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, rhs):
        return ModelDB(list.__add__(self, rhs))

    def __delitem__(self, i):
        del self.models[i]
        self.remove_unused()

    def __contains__(self, item):
        return item in self.models

    def __setitem__(self, i, value):
        self.models[i] = value

    def append(self, model):
        self.models.append(model)

    def append_lst(self, lst, constraints=None):
        self.models += lst
        if constraints:
            self.constrain(constraints)

    def get_all_abs(self, iden, param, locked=False, constraints=None):
        """return a list of desired param values from all models"""
        x = []
        y = []
        if constraints is not None:
            self.constrain(constraints)
        lst = self.models
        if locked:
            for item in lst:
                abslist = item.absorber_list
                ab = abslist.get_item(iden)
                if ab.locked(param):
                    x.append(float(getattr(ab, param)))
                    y.append(float(item.chi2))
        else:
            for item in lst:
                abslist = item.absorber_list
                ab = abslist.get_item(iden)
                if ab is None:
                    raise Exception(item.xmlfile + " has no absorber: " + iden)
                try:
                    x.append(float(getattr(ab, param)))
                    y.append(float(item.chi2))
                except:
                    pass

        if len(x) == 0 or len(y) == 0 or len(x) != len(y):
            raise Exception("ill condittioned input: \n  x=%s\n  y=%s" % (str(x), str(y)))

        return x, y

    def get_attr_lst(self, attr, cond_fn, *args):
        """ 
        cond_fn should be some boolean function definition
        """
        if args:  # if attr is actually a function
            return [getattr(it, attr)(*args) for it in self.models if cond_fn(it)]
        else:
            return [getattr(it, attr) for it in self.models if cond_fn(it)]

    def append_db(self, dbfile):
        """appends another xmldb from filename to the current db"""
        tree = et.parse(dbfile)
        root = tree.getroot()

        for key in Model.model_classes.keys():
            parent = root.find(key + 's')
            objlist = data_types.ObjList.list_from_xml(parent)  # instantiate all absorber/contpoint/view/etc
        model_list = []
        models = root.find('ModelDB').findall('model')

        for model in models:
            # get model data (includes an id mapping to something in objlisr)
            kwargs = {}
            for key, val in dict(model.attrib).items():
                kwargs[key] = val
            try:
                model_list.append(Model(**kwargs))
            except:
                raise Exception(str(kwargs))

    @staticmethod
    def filter(inp, filters):
        for key, val in filters.items():
            if key == 'chi2':
                inp = [item for item in inp if item.chi2 < val]
            else:
                if type(val) is dict:
                    for param, rng in val.items():
                        if param == 'shift':
                            inp = [item for item in inp if rng[0] <= item.get_shift(key, 'H') <= rng[-1]]
                        elif type(rng) in [str, bool]:
                            inp = [item for item in inp if item.get_datum(key, 'Absorber', param) == rng]
                        else:
                            inp = [item for item in inp if
                                   rng[0] <= float(item.get_datum(key, 'Absorber', param)) <= rng[-1]]
                else:
                    inp = [item for item in inp if getattr(item, key) == val]
        return inp

    def constrain(self, constraints):
        """
        example constraints:   
            constraints={"chi2":123,"params":3,"pixels":2345,"D":{"N":(12.3,14.3),"b":(15,16)}}
        """

        if type(constraints) is dict:
            constraints = Constraint(**constraints)

        for item in [it for it in self.models if not it in constraints]:
            self.remove(item)

    def get(self, xmlfile, chi2, pixels, params=None):
        """get from xml fit file"""
        mod = Model(xmlfile=xmlfile, chi2=chi2, pixels=pixels, params=params)
        # test for unphysical values
        mod.check_vals()
        self.models.append(mod)

    def get_lst_from_id(self, iden, attr):
        """get all models with a given continuum"""
        return [item for item in self.models if getattr(item, attr) == iden]

    def get_best_lst(self, iden=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.models = sorted(self.models, key=lambda x: x.chi2)
            return [(mod, mod.chi2) for mod in self.models]
        else:
            return self.get_locked(iden, param)  # already sorted

    def get_locked(self, iden, tag, param):
        tmp = []

        for mod in self.models:
            try:
                if to_bool(mod.get_datum(iden, tag, param + "Locked")):
                    tmp.append(mod)
            except:
                pass
        tmp = sorted(tmp, key=lambda x: x.chi2)
        return [(mod.get_datum(iden, tag, param), mod.chi2) for mod in tmp]

    def get_vel(self, id1, id2):
        return [item.get_vel(id1, id2) for item in self.models]

    def grab(self, xmlfile, chi2, pixels, params, **kwargs):
        """grab from xml file"""
        # need to re-instantiate xml file
        mod = Model(xmlfile=xmlfile, chi2=chi2, pixels=pixels, params=params, **kwargs)
        self.models.append(mod)
        return

    @staticmethod
    def read(filename, returndb=True, verbose=False, print_time=False):
        """read from xml db, return inputs for Model"""

        if not filename.endswith('.xml'):
            try:
                return ModelDB.load_models(filename)
            except:
                msg = "ModelDB.read: input file %s is" % (filename)
                msg += " of incorrect format. Please verify"
                raise Exception(msg)

        import time
        t = time.time()
        tree = et.parse(filename)
        t0 = time.time()
        if print_time:
            print("time to parse file: %lf" % (t0 - t))
        root = tree.getroot()
        # root=xmlutils.Model_xml.get_root(filename)
        # build all data first before instantiating individual models
        for key in list(Model.model_classes.keys()):
            if verbose:
                print('reading ' + key + 's instances...')
            parent = root.find(key + 's')
            objlist = data_types.ObjList.list_from_xml(parent,
                                                       verbose)
            # instantiate all absorber/contpoint/view/etc data.

        t1 = time.time()
        if time:
            print("time to get objlist: %lf" % (t1 - t0))

        model_list = []
        models = root.find('ModelDB').findall('model')

        for model in models:
            # get model data (includes an id mapping to something in objlisr)
            kwargs = {}
            for key, val in dict(model.attrib).items():
                kwargs[key] = val
            try:
                model_list.append(Model(**kwargs))
            except:
                raise Exception(str(kwargs))

        t2 = time.time()
        if time: print("instantiate models: %lf" % (t2 - t1))

        if len(model_list) == 0:
            raise Exception("no models saved")

        if not returndb:
            return model_list
        else:
            return ModelDB(models=model_list)

    @staticmethod
    def merge(db1, db2):
        if not type(db1) is type(db2):
            raise TypeError("conflicting types: got %s and %s" % (str(type(db1)), str(type(db2))))
        if type(db1) is str:
            db1, db2 = ModelDB.load_models(db1), ModelDB.load_models(db2)
        elif not type(db1) is ModelDB:
            msg = "expected either string or ModelDB instance.  instead got type %s" % (str(type(db1)))
            raise TypeError(msg)
        else:
            pass

        db1.models += db2.models

    def set_val(self, model, **kwargs):
        """set the values of a given data element"""
        if type(model) is Model:
            for item in self.models:
                if item is model:
                    item.set_val(**kwargs)
                    return
            raise Exception("model.ModelDB.set_val:  model not found")
        raise Exception('model should be Model type, not string')
        # elif type(model) is str:
        #     model = self.get_model(iden)
        #     model.set_val(**kwargs)

    @staticmethod
    def dump_models(db, fname=None):
        if not fname:
            fname = db.name
        if not fname.endswith(".obj"):
            fname += ".obj"
        with open(fname, "wb") as f:
            pickle.dump(db, f)

    @staticmethod
    def load_models(fname):
        with open(fname, "rb") as f:
            db = pickle.load(f)
        return db


# some helper functions:
# ------------------------------

def inv_dict(tag, dic=Model.model_classes):
    if tag in dic.values():
        tmp = {v: k for k, v in dic.items()}
        return tmp[tag]


def to_bool(string):
    string = string.lower()
    if string == "true":
        return True
    else:
        return False


def check_for_conflicts(root):
    """for some unknown reason, duplicates of a certain node will be printed.  
    this does not fix the underlying cause, but is a fix to prevent duplicate 
    printing.  If two nodes have the same id, but differing contents, an 
    exception will be raised"""
    ids = []
    for item in root:
        try:
            iden = item.get("id")
            assert (iden not in ids)
        except AssertionError:
            # code
            pass


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    copied and modified from:
    http://stackoverflow.com/questions/17402323/use-xml-etree-elementtree-to-write-out-nicely-formatted-xml-files
    """
    from xml.dom import minidom
    rough_string = et.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
