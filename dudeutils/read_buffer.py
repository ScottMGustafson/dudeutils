import model
import io
import subprocess
        
def populate_database(abs_ids=[],keep=False,
        separator="\n\n", path=None, db=None, constraints=None, buff=None):
    """
    parse all xmlfiles in a given directory and returns a ModelDB instance

    Input:
    ------
    path : specified path to check

    Output:
    -------
    ModelDB instance

    Raises:
    -------
    Exception

    """
    if buff:  #passes buffer to model.Model.read(), then to data_types.read(), 
              #then to xml.etree.ElementTree.parse() as io.BytesIO


        strlst=buff.decode().strip(separator).split(separator)
        buff=[io.BytesIO(item.replace("\n","").encode()) for item in strlst]
        models = [model.Model(buff=item) for item in buff]

    else:
        if type(abs_ids)!=list:
            raise Exception("abs_ids needs to be list")

        if len(abs_ids)==0:
            raise Exception("need at least one absorber specified.  abs_ids=[]")

        #read files on disk
        files=[]
        if not path:
            path=os.getcwd()
        for f in os.listdir(path):
            if split(f)[-1].startswith('iteration') and split(f)[-1].endswith(".xml"):
                files.append(f)

        if len(files)==0 and not db:
            raise Exception("too few files")

        models=[Model( xmlfile=join(path,f), abs_ids=abs_ids ) for f in files]
        if not keep:
            for f in files:
                os.remove(join(path,f))

    if bool(db): #if None or []
        if type(db) is str:
            db=load_from_db(db)
        db.append_lst(models,constraints=constraints)
    else:
        db=model.ModelDB(models=models,constraints=constraints)


    return db


def run_optimize(fname,step=False, verbose=False, 
        to_buffer=True, method='dude',timeout=None):
    """
    call dude and run OptimizeXML from the command line.

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
    popt, pcov:   if method!= 'dude', outputs optimized parameters and 
                  covariance matrix (tuple)

    Raises:
    -------
    None
    """
    if method=='dude':
        #call dude from the command line and call its Levenberg-marquardt algo.
        commands=["java","-cp",
                  "/home/scott/programming/dude/jd/build", 
                  "dude.commandline.OptimizeXML", 
                  fname
                  ]
        if step:
            commands.append('step')
        if to_buffer:
            commands.append('to_buffer')
        if verbose:
            print("running: %s"%(" ".join(commands)))


        #now run the optimizer
        if to_buffer:
            return subprocess.check_output(commands,timeout=timeout)
                   
        else:
            subprocess.call(commands,timeout=timeout)
            return None
    else: #use code from this project
        raise Exception("not yet fully implemented as of 2016-02-29")


if __name__ == "__main__":
    buff=run_optimize("/home/scott/research/J0744+2059/test.xml",step=True, to_buffer=True)
    db=populate_database(buff=buff)
    print(len(db))

