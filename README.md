 # Dudeutils

A framework for dude spectrum models and .fits data

## Important data structures
`dudeutils.data_types.Data`:  Parent class for `Absorber`, `ContinuumPoint`, `SingleView`, `VelocityView` and `Region` classes.

`dudeutils.data_types.ObjList`: container for `data_types.Data` subclasses.  Elements are stored in `ObjList._pool`, which you as a user should never directly use, unless absolutely necessary.

`dudeutils.model.Model`:  class for model data.  includes functions for parsing, dumping and manipulating dude models.

`dudeutils.model.ModelDB`: container class for Model.  Acts like an extended list.

`dudeutils.spec_parser.Spectrum`: parent class for `TextSpectrum` and `FitsSpectrum`.  Provides utilities for interfacing with data and fitting model data.   


## example usage
For more detailed usage, see the 'Software' chapter in my PhD Thesis

```python
from dudeutils.model import Model, ModelDB
```


to get a single model from `db`
```python
model=db.get_model(id='the_model')
```

reset all data on a model to another one:
```python
model.get_model('another_file.xml')
```
get an individual absorber/continuum point/other datum
```python
datum = model.pop('the_id','Absorber')
```

append datum:
```python
model.append(datum)
```

set a new value for a datum:
```python
new_vals = {'N':21.5,'b':17.2,'z':1.23,'bLocked':True,'zLocked':False,'ionName':'H I','id':'new_id'}
model.set('old_id','Absorber',**new_vals)
```

## Random sampling of models
run from random_sampling module:
```python
from dudeutils.random_sampling import run
from dudeutils.model_csv_io import ModelIO

run('/path/to/config.cfg')
```

If you need to tie two absorbers' parameters (i.e. redshift), use the `tie` parameter:

```python
run('/path/to/config.cfg', tie=('H_z', 'D_z'))
```
Models will be written as a csv with headers following the format:  {absorber id}_{absorber parameter}, i.e. the redshift of an absorber named "abc" will be `abc_z`

#### be sure to write all of your changes back to file if you want them to take effect:

write the model back to the original xml

```python
Model.write()
```

Saving an entire db:

```python
db.write()
```

The models are also pickleable, which is often more useful.  We have built-in functions for this too.

```python
ModelDB.dump_models(db,'name.obj')
db=ModelDB.load_models('name.obj')
```

#### XML formats
A large portion of this software boils down to parsing and storing xml data according to a few specific formats using xml.etree.

##### fit database formats:
each distinct data type will be stored separately with a unique id as an identifier.

for example:
```xml
<AbsorberList id="abslistID">
    <Absorber N="13.4" NError="0.0" NLocked="true" b="11.34" bError="0.0" bLocked="true" id="anAbsorber" ionName="C III" z="2.9" zError="0.0" zLocked="true" />
</AbsorberList>
<ContinuumPointList id="contpntID">
    <ContinuumPoint id="pnt" x="1234.5" xError="0.0" xLocked="true" y="4.0E-14" yError="0.0" yLocked="true" />
</ContinuumPointList>
<ModelDB>
    <Model id="modelid", AbsorberList="abslistID", ContinuumPointList="contpntID", chi2="1892", pixels="187", params="12">
</ModelDB>
```

##### an individual fit's xml file
This should be the standard dude xml output.  
    
    
