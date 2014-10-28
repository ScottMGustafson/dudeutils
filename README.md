#Dudeutils

A few scripts to assist in data analysis with our group's in-house spectra fitting software.
Basically it parses and edits the xml file dude reads from to display and interact with your spectrum.  Not very elegant as it requires lots of saving and reloading on the Dude end of things, but I didn't want to mess around with a piece of software that already works well.

##example usage:

###Whilst running dude:
 1. start a python3 session in the terminal
 2. make and save whatever changes to your spectrum on dude.
 3. perform your action in the terminal
 4. if you made any changes to the spectrum and want to see them take effect, reload the xml file in dude
 5. repeat steps 2-4 as necessary
 

####start a new Model_db from a given xml fit, which will give the first model:
```python
import data_structures as data
chi2=123
pixels=1234
dbfile=None #if None, one will be created for you
db = data.newdb('xmlfile.xml',chi2,pixels,dbfile=None)
```

####access an alread-existing Model_db from file:
returns a ModelDB instance
```python
db = data.ModelDB.read('the_file.xml')
```

or to get as a list of models:
```python
model_list = data.ModelDB.read('the_file.xml',return_db=False) 
```

####manipulations on a single model:
    

to get a single model from `db`
```python
model=db.get_model(id='the_model')
```
or to grab and append from another xml fit:
```python
db.grab(id='the_model')
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

randomly set a value: this example randomly sets the x value between 3403 and 3405 angstroms
```python
data.random('xmlfile.xml','cont_pt_id','ContinuumPoint','x',[3403.,3405.],"id_for_model")
```

####be sure to write all of your changes back to file if you want them to take effect:

write the model back to the original xml

```python
model.write()
```

Saving an entire db:

```python
db.write()
```

####XML formats
A large portion of this software boils down to parsing and storing xml data according to a few specific formats using xml.etree.

#####fit database formats:
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

#####an individual fit's xml file
This should be the standard dude xml output.  If you have no idea whe this means, then you don't have dude and this entire repository is probably not for you
    
    
