import os
import subprocess
import xml.etree.ElementTree as et
from os.path import split, join

from dudeutils.model import Model, ModelDB


def run_optimize(fname, step=False, verbose=False,
                 method='dude', timeout=None, **kwargs):
    """
    call dude and run commandline.OptimizeXML from the command line.

    input:
    ------
    fname:  source dude xml file (string)
    step (bool): if True, calls dude's step-iteration procedure,
        dumping some number of fits for each step of the procedure, named 
        \"iteration_%d.xml\"

    ouput:
    ------
    popt, pcov:   if method != 'dude', outputs optimized parameters and 
                  covariance matrix (tuple)

    Raises:
    -------
    None
    """
    if timeout:
        timeout = float(timeout)
    if method == 'dude':
        run_dude_optimizer(fname, step=step, verbose=verbose, timeout=timeout)
    else:  # use code from this project
        raise NotImplementedError("doesn't yet work as of 2016-02-29")
        # model = Model(xmlfile=fname)
        # src_data = model.flux
        # popt, pcov = optimizer.optimize(src_data, model)
        return popt, pcov

def run_dude_optimizer(fname, step=False, verbose=False, timeout=None):
    # call dude from the command line and call its Levenberg-marquardt algo.
    commands = ["java", "-cp",
                "/home/scott/programming/dude/jd/build",
                "dude.commandline.OptimizeXML",
                fname
                ]
    if step:
        commands.append('step')
    if verbose:
        print("running: %s" % (" ".join(commands)))
    # now run the optimizer
    return subprocess.check_output(commands, timeout=timeout)

def run_sim_annealing_optimizer(fname, step=False, verbose=False, timeout=None):
    # call dude from the command line and call its Levenberg-marquardt algo.

    return

def newdb(xmlfile, dbfile=None, params=None, **kwargs):
    """get a model, append to new database"""
    model = Model(xmlfile=xmlfile)
    return ModelDB(dbfile, [model], **kwargs)


def get_model(xmlfile, chi2=0., pixels=0.):
    """
    get a single model instance from an xmlfile.  This exists just to provide 
    multiple points of entry to get a model

    Input:
    ------
    xmlfile : a dude xml fit file
    chi2 : chi-square.  default=0.
    pixels : number of pixels being optimized.  default=0

    Output:d
    -------
    model.Model instance 

    """
    return Model(xmlfile=xmlfile, chi2=chi2, pixels=pixels)


def populate_database(abs_ids=None, separator="\n\n", return_list=False, xml_lst=None,
                      path=None, db=None, constraints=None):
    """
    parse all xmlfiles in a given directory and returns a ModelDB instance

    Input:
    ------
    abs_ids: list of absorbers ids to keep (str)
    separator: string separator to separate models
    path : specified path to check if reading from file
    db: if appending to specified database.
    constraints:  model constraints.  (see constraints.Constraints)
    xml_lst: list of xml strings to parse rather than reading from disk

    Output:
    -------
    ModelDB instance

    Raises:
    -------
    Exception

    """


    if isinstance(xml_lst, list):
        models = [Model(xmlfile=xml, abs_ids=abs_ids) for xml in xml_lst]

    else:
        # read files on disk
        files = []
        if not path:
            path = os.getcwd()
        for f in os.listdir(path):
            if split(f)[-1].startswith('iteration') and split(f)[-1].endswith(".xml"):
                files.append(f)

        if len(files) == 0 and not db:
            raise Exception("too few files")
        models = []
        for f in files:
            try:
                models.append(Model(xmlfile=join(path, f), abs_ids=abs_ids))
            except:
                print(join(path, f), " failed to parse")
        for f in files:
            os.remove(join(path, f))
    if return_list:
        return models
    elif db:  # if None or []
        if type(db) is str:
            db = load_from_db(db)
        db.append_lst(models, constraints=constraints)
        return db
    else:
        return ModelDB(models=models, constraints=constraints)


def append_db(db, abs_ids, keep=False, path=None, constraints=None):
    """populate existing database"""
    populate_database(abs_ids, path=path, db=db, constraints=constraints)


def dump_models(db, fname=None):
    """alias to ModelDB.dump_models"""
    ModelDB.dump_models(db, fname)


def load_models(fname):
    """alias to ModelDB.load_models"""
    return ModelDB.load_models(fname)


def load_from_db(dbxml):
    """alias to load an xml database"""
    return ModelDB.read(dbxml)


def get_db(dbfile):
    """alias to load an xml database"""
    return ModelDB.read(dbfile)


def parse(filename, ab_ids, path='/home/scott/research/J0744+2059/'):
    """clear a db file of all absorbers not in ab_ids"""

    tree = et.parse(join(path, filename))

    root = tree.getroot()
    parent = root.find('AbsorberLists')
    models = root.find('ModelDB').findall('model')

    for ablist in parent.findall('AbsorberList'):
        for ab in ablist.findall('Absorber'):
            if ab.get('id') not in ab_ids:
                ablist.remove_item(ab)
        if len(list(ablist)) == 0:
            parent.remove_item(ablist)
    tree.write(join(path, filename))
