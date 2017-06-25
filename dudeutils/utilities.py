from dudeutils.model import Model, ModelDB
import os
from os.path import split, join
import subprocess
import xml.etree.ElementTree as et
import io


def run_optimize(fname, step=False, verbose=False, to_buffer=False,
                 method='dude', timeout=None, **kwargs):
    """
    call dude and run commandline.OptimizeXML from the command line.

    input:
    ------
    fname:  source dude xml file (string)
    step (bool): if True, calls dude's step-iteration procedure,
        dumping some number of fits for each step of the procedure, named 
        \"iteration_%d.xml\"
    to_buffer (bool): if True, writes dude output to buffer instead of disk.  
        models are then passed to xml.etree.ElementTree as BytesIO or StringIO

    ouput:
    ------
    popt, pcov:   if method != 'dude', outputs optimized parameters and 
                  covariance matrix (tuple)

    Raises:
    -------
    None
    """
    if timeout: timeout = float(timeout)
    if method == 'dude':
        # call dude from the command line and call its Levenberg-marquardt algo.
        commands = ["java", "-cp",
                    "/home/scott/programming/dude/jd/build",
                    "dude.commandline.OptimizeXML",
                    fname
                    ]
        if step:
            commands.append('step')
        if to_buffer:
            commands.append('to_buffer')
        if verbose:
            print("running: %s" % (" ".join(commands)))

        # now run the optimizer
        return subprocess.check_output(commands, timeout=timeout)
    else:  # use code from this project
        raise Exception("doesn't yet work as of 2016-02-29")
        model = Model(xmlfile=fname)
        src_data = model.flux
        popt, pcov = optimizer.optimize(src_data, model)
        return popt, pcov


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


def populate_database(abs_ids=None, separator="\n\n", return_list=False,
                      path=None, db=None, constraints=None, buff=None):
    """
    parse all xmlfiles in a given directory and returns a ModelDB instance

    Input:
    ------
    abs_ids: list of absorbers ids to keep (str)
    separator: string separator to separate models in buffer stream
    path : specified path to check if reading from file
    db: if appending to specified database.
    constraints:  model constraints.  (see constraints.Constraints)
    buff: a buffer to read in models from memory

    Output:
    -------
    ModelDB instance

    Raises:
    -------
    Exception

    """
    if buff:  # passes buffer to Model.read(), then to data_types.read(),
        # then to xml.etree.ElementTree.parse() as io.BytesIO

        strlst = buff.decode().strip(separator).split(separator)
        buff = [io.BytesIO(item.replace("\n", "").encode()) for item in strlst]
        models = [Model(buff=item) for item in buff]

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
    populate_database(abs_ids, keep=keep, path=path, db=db, constraints=constraints)


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
                ablist.remove(ab)
        if len(list(ablist)) == 0:
            parent.remove(ablist)
    tree.write(join(path, filename))


if __name__ == '__main__':
    pass
